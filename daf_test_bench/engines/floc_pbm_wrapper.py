"""
Floc Kinetics PBM Engine Wrapper

Facade/Adapter for the Floc Kinetics Population Balance Model.
Implements the IDAFPlant interface for integration with the Test Harness.

This wrapper translates standardized test configurations into PBM-specific
API calls and extracts results in the standardized BenchmarkResult format.
"""

import time
import traceback
import logging
from typing import Dict, Any, Optional, List
import numpy as np

# Import PBM dependencies
try:
    from floc_kinetics_pbm.src.pbe_solver import PopulationBalanceModel
    from floc_kinetics_pbm.src.properties import FlocProperties
    PBM_AVAILABLE = True
except ImportError as e:
    # Fallback: mark wrapper as unavailable
    logging.warning(f"Floc PBM dependencies not available: {e}")
    PopulationBalanceModel = None
    FlocProperties = None
    PBM_AVAILABLE = False

from ..interfaces.idaf_plant import IDAFPlant
from ..data_models import BenchmarkResult, ScientificMetrics, ComputationalMetrics


class FlocPBMWrapper(IDAFPlant):
    """
    Wrapper for Floc Kinetics PBM engine.

    The Floc Kinetics PBM solves the discretized population balance equation
    for floc aggregation and breakage in 0D (well-mixed) or as a local operator
    in CFD simulations.
    """

    def __init__(self):
        """Initialize Floc Kinetics PBM wrapper."""
        super().__init__(engine_name="Floc_Kinetics_PBM")
        self._pbm_model: Optional['PopulationBalanceModel'] = None
        self._properties: Optional[FlocProperties] = None
        self._initial_psd: Optional[np.ndarray] = None
        self._time_points: Optional[np.ndarray] = None
        self._psd_history: Optional[np.ndarray] = None
        self._final_psd: Optional[np.ndarray] = None
        self._statistics: Optional[Dict[str, Any]] = None
        self._start_time: float = 0.0
        self._end_time: float = 0.0
        self._simulation_time: float = 0.0
        self._num_timesteps: int = 0
        self._shear_rate: float = 50.0

        # Set up logging
        self.logger = logging.getLogger(f"{__name__}.{self.engine_name}")

    def setup(self, configuration: Dict[str, Any]) -> bool:
        """
        Configure the PBM simulation.

        Expected configuration structure:
            stimuli:
                pbm:
                    num_size_classes: [int] Number of bins
                    d_min: [m] Minimum diameter
                    d_max: [m] Maximum diameter
                    initial_mean_diameter: [m]
                    initial_concentration: [#/m³]
                    kernel_coeffs:
                        collision_efficiency: [-]
                        shear_rate: [1/s]
                        temperature: [K]
                        primary_diameter: [m]
                        fractal_dimension: [-]
                        binding_strength: [N/m²]
                    simulation:
                        total_time: [s]
                        num_timesteps: [-]

        Args:
            configuration: Test configuration dictionary

        Returns:
            bool: True if setup successful
        """
        try:
            if not PBM_AVAILABLE:
                self.logger.error("PBM dependencies not available")
                return False

            self._config = configuration
            self.logger.info(f"Setting up {self.engine_name}")
            self.logger.info(f"Test: {configuration.get('test_name', 'Unknown')}")

            # Extract PBM stimuli
            stimuli = configuration.get('stimuli', {})
            pbm_params = stimuli.get('pbm', {})

            if not pbm_params:
                # Use default parameters
                pbm_params = self._get_default_parameters()
                self.logger.info("Using default PBM parameters")

            # Size parameters
            num_classes = pbm_params.get('num_size_classes', 20)
            d_min = pbm_params.get('d_min', 1e-6)
            d_max = pbm_params.get('d_max', 1e-3)

            # Kernel coefficients
            kernel_coeffs = pbm_params.get('kernel_coeffs', {})
            self._shear_rate = kernel_coeffs.get('shear_rate', 50.0)

            self._properties = FlocProperties(
                collision_efficiency=kernel_coeffs.get('collision_efficiency', 0.1),
                temperature=kernel_coeffs.get('temperature', 293.15),
                primary_particle_diameter=kernel_coeffs.get('primary_diameter', 1e-6),
                fractal_dimension=kernel_coeffs.get('fractal_dimension', 2.3),
                binding_strength=kernel_coeffs.get('binding_strength', 1e3)
            )

            # Create PBM model
            self._pbm_model = PopulationBalanceModel(
                n_bins=num_classes,
                d_min=d_min,
                d_max=d_max,
                properties=self._properties
            )

            # Initial PSD
            if 'influent_psd' in pbm_params:
                self._initial_psd = np.array(pbm_params['influent_psd'])
                if len(self._initial_psd) != self._pbm_model.n_bins:
                    raise ValueError(
                        f"PSD length ({len(self._initial_psd)}) must match "
                        f"size classes ({self._pbm_model.n_bins})"
                    )
            else:
                # Default: Log-normal distribution centered at specified diameter
                mean_size = pbm_params.get('initial_mean_diameter', 10e-6)
                std_size = pbm_params.get('initial_std_diameter', 3e-6)
                total_conc = pbm_params.get('initial_concentration', 1e13)
                self._initial_psd = self._create_lognormal_psd(
                    self._pbm_model.diameters, mean_size, std_size, total_conc
                )

            # Simulation parameters
            sim_params = pbm_params.get('simulation', {})
            self._simulation_time = sim_params.get('total_time', 600.0)  # 10 min default
            self._num_timesteps = sim_params.get('num_timesteps', 100)

            self.logger.info(f"Size classes: {self._pbm_model.n_bins}")
            self.logger.info(f"Size range: {self._pbm_model.d_min*1e6:.2f} - {self._pbm_model.d_max*1e6:.0f} μm")
            self.logger.info(f"Initial total concentration: {np.sum(self._initial_psd):.2e} #/m³")
            self.logger.info(f"Collision efficiency: {self._properties.collision_efficiency}")
            self.logger.info(f"Shear rate: {self._shear_rate:.1f} 1/s")
            self.logger.info(f"Simulation time: {self._simulation_time:.0f} s")

            self._is_setup = True
            return True

        except Exception as e:
            self.logger.error(f"ERROR in setup: {str(e)}")
            self.logger.debug(traceback.format_exc())
            return False

    def initialize(self) -> bool:
        """
        Initialize the PBM simulation.

        Precomputes aggregation and breakage kernel matrices.

        Returns:
            bool: True if initialization successful
        """
        try:
            self.validate_workflow("initialize")
            self.logger.info("Initializing PBM solver")

            # Verify configuration
            if self._pbm_model is None or self._initial_psd is None:
                raise RuntimeError("PBM model and initial PSD must be configured")

            # Log mass conservation check
            initial_mass = np.sum(self._initial_psd * (self._pbm_model.volumes))
            self.logger.info(f"Initial total volume: {initial_mass:.2e} m³/m³")

            self._is_initialized = True
            self.logger.info("Initialization complete")
            return True

        except Exception as e:
            self.logger.error(f"ERROR in initialize: {str(e)}")
            self.logger.debug(traceback.format_exc())
            return False

    def run(self) -> bool:
        """
        Run the PBM simulation.

        Solves the transient PBE from t=0 to t=total_time.

        Returns:
            bool: True if simulation successful
        """
        try:
            self.validate_workflow("run")
            self.logger.info("Running PBM simulation")
            self.logger.info(f"Integrating from t=0 to t={self._simulation_time:.0f} s")

            self._start_time = time.time()

            # Time points for integration
            self._time_points = np.linspace(0, self._simulation_time, self._num_timesteps + 1)

            # Solve PBE using the PopulationBalanceModel.solve() method
            self.logger.info("Solving population balance equation...")
            t, n_solution = self._pbm_model.solve(
                n_initial=self._initial_psd,
                t_span=self._time_points,
                G=self._shear_rate,
                recompute_kernels=True
            )

            # Store results
            self._psd_history = n_solution
            self._final_psd = n_solution[-1, :]

            # Calculate statistics using the model's summary_statistics method
            self._statistics = self._pbm_model.summary_statistics(self._final_psd)

            self._end_time = time.time()
            elapsed = self._end_time - self._start_time

            # Mass conservation check
            initial_mass = np.sum(self._initial_psd * self._pbm_model.volumes)
            final_mass = np.sum(self._final_psd * self._pbm_model.volumes)
            mass_error = abs(final_mass - initial_mass) / initial_mass * 100 if initial_mass > 0 else 0.0

            self.logger.info(f"Simulation complete in {elapsed:.3f} s")
            self.logger.info(f"Final mean diameter (d50): {self._statistics['d50']*1e6:.2f} μm")
            self.logger.info(f"Final total concentration: {self._statistics['total_number']:.2e} #/m³")
            self.logger.info(f"Mass conservation error: {mass_error:.4f}%")

            self._run_completed = True
            return True

        except Exception as e:
            self.logger.error(f"ERROR in run: {str(e)}")
            self.logger.debug(traceback.format_exc())
            return False

    def get_metrics(self) -> Dict[str, Any]:
        """
        Extract metrics in standardized format.

        Returns:
            dict: BenchmarkResult dictionary
        """
        try:
            self.validate_workflow("get_metrics")

            # Determine run status
            run_status = "Success" if self._run_completed else "Failed"

            # Scientific metrics
            sci_metrics = ScientificMetrics()
            if self._statistics:
                sci_metrics.mean_floc_size_d50 = self._statistics['d50']

                # Calculate mass conservation error
                if self._initial_psd is not None and self._final_psd is not None:
                    initial_mass = np.sum(self._initial_psd * self._pbm_model.volumes)
                    final_mass = np.sum(self._final_psd * self._pbm_model.volumes)
                    mass_error = abs(final_mass - initial_mass) / initial_mass * 100 if initial_mass > 0 else 0.0
                    sci_metrics.mass_conservation_error_pct = mass_error

                # Additional PBM metrics
                sci_metrics.additional_metrics = {
                    'mean_diameter': self._statistics['mean_diameter'],
                    'd10_number': self._statistics['d10_number'],
                    'd50_number': self._statistics['d50_number'],
                    'd90_number': self._statistics['d90_number'],
                    'd10_volume': self._statistics['d10_volume'],
                    'd50_volume': self._statistics['d50_volume'],
                    'd90_volume': self._statistics['d90_volume'],
                    'total_concentration': self._statistics['total_number'],
                    'total_volume_fraction': self._statistics['volume_fraction']
                }

            # Computational metrics with real memory measurement
            try:
                import psutil
                import os
                process = psutil.Process(os.getpid())
                peak_ram_gb = process.memory_info().rss / (1024 ** 3)
            except ImportError:
                # Fallback if psutil not available
                peak_ram_gb = None

            comp_metrics = ComputationalMetrics(
                wall_clock_sec=self._end_time - self._start_time if self._start_time > 0 else None,
                cpu_hours=(self._end_time - self._start_time) / 3600.0 if self._start_time > 0 else None,
                peak_ram_gb=peak_ram_gb,
                converged=True,
                num_iterations=self._num_timesteps
            )

            # Create result
            result = BenchmarkResult.create(
                engine_name=self.engine_name,
                test_name=self._config.get('test_name', 'Unknown'),
                run_status=run_status,
                sci_metrics=sci_metrics,
                comp_metrics=comp_metrics,
                run_log=self._get_log_contents()
            )

            return result.to_dict()

        except Exception as e:
            error_log = self._get_log_contents() + f"\nERROR in get_metrics: {str(e)}\n{traceback.format_exc()}"
            result = BenchmarkResult.create(
                engine_name=self.engine_name,
                test_name=self._config.get('test_name', 'Unknown') if self._config else 'Unknown',
                run_status="Failed",
                run_log=error_log
            )
            return result.to_dict()

    def get_field_data(self, variable_name: str) -> Optional[np.ndarray]:
        """
        Get field data.

        Args:
            variable_name: Name of variable
                - "psd" or "particle_size_distribution": Final PSD
                - "psd_history": Complete time history
                - "size_classes": Size class boundaries

        Returns:
            Requested data array, or None if not available
        """
        if variable_name in ['psd', 'particle_size_distribution']:
            return self._final_psd.copy() if self._final_psd is not None else None
        elif variable_name == 'psd_history':
            return self._psd_history.copy() if self._psd_history is not None else None
        elif variable_name == 'size_classes':
            return self._pbm_model.diameters.copy() if self._pbm_model is not None else None
        elif variable_name == 'time_points':
            return self._time_points.copy() if self._time_points is not None else None

        return None

    def finalize(self) -> None:
        """Clean up resources."""
        self.logger.info("Finalizing PBM engine")
        # Clear large arrays to free memory
        self._psd_history = None

    def _create_lognormal_psd(
        self,
        size_classes: np.ndarray,
        mean: float,
        std: float,
        total_conc: float
    ) -> np.ndarray:
        """
        Create a log-normal particle size distribution.

        Args:
            size_classes: Size class diameters [m]
            mean: Mean diameter [m]
            std: Standard deviation [m]
            total_conc: Total particle concentration [#/m³]

        Returns:
            PSD array [#/m³]
        """
        # Correct log-normal conversion
        log_mean = np.log(mean)
        log_std = np.sqrt(np.log(1 + (std / mean)**2)) if mean > 0 else 0

        psd = np.exp(-0.5 * ((np.log(size_classes) - log_mean) / log_std) ** 2)
        psd = psd / np.sum(psd) * total_conc

        return psd

    def _get_log_contents(self) -> str:
        """
        Get logging contents.

        In the future, this could return contents from a logging handler.
        For now, returns a simple summary.
        """
        return f"PBM simulation log for {self.engine_name}"

    def _get_default_parameters(self) -> Dict[str, Any]:
        """Get default PBM parameters."""
        return {
            'num_size_classes': 20,
            'd_min': 1e-6,
            'd_max': 1e-3,
            'initial_mean_diameter': 10e-6,
            'initial_std_diameter': 3e-6,
            'initial_concentration': 1e13,
            'kernel_coeffs': {
                'collision_efficiency': 0.1,
                'shear_rate': 50.0,
                'temperature': 293.15,
                'primary_diameter': 1e-6,
                'fractal_dimension': 2.3,
                'binding_strength': 1e3
            },
            'simulation': {
                'total_time': 600.0,
                'num_timesteps': 100
            }
        }
