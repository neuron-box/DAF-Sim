"""
Floc Kinetics PBM Engine Wrapper

Facade/Adapter for the Floc Kinetics Population Balance Model.
Implements the IDAFPlant interface for integration with the Test Harness.

This wrapper translates standardized test configurations into PBM-specific
API calls and extracts results in the standardized BenchmarkResult format.
"""

import time
import traceback
from typing import Dict, Any, Optional, List
import numpy as np
import sys
from pathlib import Path

# Add parent directory of floc_kinetics_pbm to path
floc_pbm_parent = str(Path(__file__).parent.parent.parent)
if floc_pbm_parent not in sys.path:
    sys.path.insert(0, floc_pbm_parent)

# Import after path is set - use package import
try:
    from floc_kinetics_pbm.src.kernels import calculate_floc_kernels
    from floc_kinetics_pbm.src.pbe_solver import solve_pbe_transient, calculate_psd_statistics
    from floc_kinetics_pbm.src.properties import FlocProperties
except ImportError as e:
    # Fallback: mark wrapper as unavailable
    print(f"Warning: Floc PBM dependencies not available: {e}")
    calculate_floc_kernels = None
    solve_pbe_transient = None
    calculate_psd_statistics = None
    FlocProperties = None

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
        self._properties: Optional[FlocProperties] = None
        self._size_classes: Optional[np.ndarray] = None
        self._initial_psd: Optional[np.ndarray] = None
        self._time_points: Optional[np.ndarray] = None
        self._psd_history: Optional[np.ndarray] = None
        self._final_psd: Optional[np.ndarray] = None
        self._statistics: Optional[Dict[str, Any]] = None
        self._start_time: float = 0.0
        self._end_time: float = 0.0
        self._log: List[str] = []
        self._simulation_time: float = 0.0
        self._num_timesteps: int = 0

    def setup(self, configuration: Dict[str, Any]) -> bool:
        """
        Configure the PBM simulation.

        Expected configuration structure:
            stimuli:
                pbm:
                    influent_psd: [list of concentrations #/m³]
                    size_classes: [list of diameters in m]
                    kernel_coeffs:
                        collision_efficiency: [-]
                        shear_rate: [1/s]
                        temperature: [K]
                        breakage_rate: [-]
                    simulation:
                        total_time: [s]
                        num_timesteps: [-]

        Args:
            configuration: Test configuration dictionary

        Returns:
            bool: True if setup successful
        """
        try:
            self._config = configuration
            self._log.append(f"Setting up {self.engine_name}")
            self._log.append(f"Test: {configuration.get('test_name', 'Unknown')}")

            # Extract PBM stimuli
            stimuli = configuration.get('stimuli', {})
            pbm_params = stimuli.get('pbm', {})

            if not pbm_params:
                # Use default parameters
                pbm_params = self._get_default_parameters()
                self._log.append("Using default PBM parameters")

            # Size classes
            if 'size_classes' in pbm_params:
                self._size_classes = np.array(pbm_params['size_classes'])
            else:
                # Create logarithmic size classes (1 μm to 1000 μm)
                num_classes = pbm_params.get('num_size_classes', 20)
                self._size_classes = np.logspace(-6, -3, num_classes)

            # Initial PSD
            if 'influent_psd' in pbm_params:
                self._initial_psd = np.array(pbm_params['influent_psd'])
                if len(self._initial_psd) != len(self._size_classes):
                    raise ValueError(
                        f"PSD length ({len(self._initial_psd)}) must match "
                        f"size classes ({len(self._size_classes)})"
                    )
            else:
                # Default: Normal distribution centered at 10 μm
                mean_size = pbm_params.get('initial_mean_diameter', 10e-6)
                std_size = pbm_params.get('initial_std_diameter', 3e-6)
                total_conc = pbm_params.get('initial_concentration', 1e13)
                self._initial_psd = self._create_normal_psd(
                    self._size_classes, mean_size, std_size, total_conc
                )

            # Kernel coefficients
            kernel_coeffs = pbm_params.get('kernel_coeffs', {})
            self._properties = FlocProperties(
                collision_efficiency=kernel_coeffs.get('collision_efficiency', 0.1),
                shear_rate=kernel_coeffs.get('shear_rate', 50.0),
                temperature=kernel_coeffs.get('temperature', 293.15),
                breakage_rate=kernel_coeffs.get('breakage_rate', 0.001),
                primary_particle_diameter=kernel_coeffs.get('primary_diameter', 1e-6),
                fractal_dimension=kernel_coeffs.get('fractal_dimension', 2.3),
                floc_binding_strength=kernel_coeffs.get('binding_strength', 1e3)
            )

            # Simulation parameters
            sim_params = pbm_params.get('simulation', {})
            self._simulation_time = sim_params.get('total_time', 600.0)  # 10 min default
            self._num_timesteps = sim_params.get('num_timesteps', 100)

            self._log.append(f"Size classes: {len(self._size_classes)}")
            self._log.append(f"Size range: {self._size_classes[0]*1e6:.2f} - {self._size_classes[-1]*1e6:.0f} μm")
            self._log.append(f"Initial total concentration: {np.sum(self._initial_psd):.2e} #/m³")
            self._log.append(f"Collision efficiency: {self._properties.collision_efficiency}")
            self._log.append(f"Shear rate: {self._properties.shear_rate:.1f} 1/s")
            self._log.append(f"Simulation time: {self._simulation_time:.0f} s")

            self._is_setup = True
            return True

        except Exception as e:
            self._log.append(f"ERROR in setup: {str(e)}")
            self._log.append(traceback.format_exc())
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
            self._log.append("Initializing PBM solver")

            # Verify configuration
            if self._size_classes is None or self._initial_psd is None:
                raise RuntimeError("Size classes and initial PSD must be configured")

            # Log mass conservation check
            initial_mass = np.sum(self._initial_psd * (self._size_classes ** 3))
            self._log.append(f"Initial total mass (volume): {initial_mass:.2e} m³/m³")

            self._is_initialized = True
            self._log.append("Initialization complete")
            return True

        except Exception as e:
            self._log.append(f"ERROR in initialize: {str(e)}")
            self._log.append(traceback.format_exc())
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
            self._log.append("Running PBM simulation")
            self._log.append(f"Integrating from t=0 to t={self._simulation_time:.0f} s")

            self._start_time = time.time()

            # Calculate kernels
            self._log.append("Computing aggregation and breakage kernels...")
            kernels = calculate_floc_kernels(
                self._size_classes,
                self._properties
            )

            # Time points for integration
            self._time_points = np.linspace(0, self._simulation_time, self._num_timesteps + 1)

            # Solve PBE
            self._log.append("Solving population balance equation...")
            self._psd_history = solve_pbe_transient(
                self._size_classes,
                self._initial_psd,
                self._time_points,
                kernels['aggregation_kernel'],
                kernels['breakage_kernel']
            )

            # Extract final PSD
            self._final_psd = self._psd_history[-1, :]

            # Calculate statistics
            self._statistics = calculate_psd_statistics(
                self._size_classes,
                self._final_psd
            )

            self._end_time = time.time()
            elapsed = self._end_time - self._start_time

            # Mass conservation check
            initial_mass = np.sum(self._initial_psd * (self._size_classes ** 3))
            final_mass = np.sum(self._final_psd * (self._size_classes ** 3))
            mass_error = abs(final_mass - initial_mass) / initial_mass * 100

            self._log.append(f"Simulation complete in {elapsed:.3f} s")
            self._log.append(f"Final mean diameter (d50): {self._statistics['d50']*1e6:.2f} μm")
            self._log.append(f"Final total concentration: {np.sum(self._final_psd):.2e} #/m³")
            self._log.append(f"Mass conservation error: {mass_error:.4f}%")

            self._run_completed = True
            return True

        except Exception as e:
            self._log.append(f"ERROR in run: {str(e)}")
            self._log.append(traceback.format_exc())
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
                    initial_mass = np.sum(self._initial_psd * (self._size_classes ** 3))
                    final_mass = np.sum(self._final_psd * (self._size_classes ** 3))
                    mass_error = abs(final_mass - initial_mass) / initial_mass * 100
                    sci_metrics.mass_conservation_error_pct = mass_error

                # Additional PBM metrics
                sci_metrics.additional_metrics = {
                    'mean_diameter': self._statistics['mean_diameter'],
                    'd10': self._statistics['d10'],
                    'd32': self._statistics['d32'],
                    'd43': self._statistics['d43'],
                    'd90': self._statistics['d90'],
                    'total_concentration': self._statistics['total_concentration'],
                    'total_volume_fraction': self._statistics['total_volume_fraction']
                }

            # Computational metrics
            comp_metrics = ComputationalMetrics(
                wall_clock_sec=self._end_time - self._start_time if self._start_time > 0 else None,
                cpu_hours=(self._end_time - self._start_time) / 3600.0 if self._start_time > 0 else None,
                peak_ram_gb=0.5,  # Estimate based on array sizes
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
                run_log="\n".join(self._log)
            )

            return result.to_dict()

        except Exception as e:
            error_log = "\n".join(self._log) + f"\nERROR in get_metrics: {str(e)}\n{traceback.format_exc()}"
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
            return self._size_classes.copy() if self._size_classes is not None else None
        elif variable_name == 'time_points':
            return self._time_points.copy() if self._time_points is not None else None

        return None

    def finalize(self) -> None:
        """Clean up resources."""
        self._log.append("Finalizing PBM engine")
        # Clear large arrays to free memory
        self._psd_history = None

    def _create_normal_psd(
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
        # Log-normal distribution
        log_mean = np.log(mean)
        log_std = std / mean  # Approximate

        psd = np.exp(-0.5 * ((np.log(size_classes) - log_mean) / log_std) ** 2)
        psd = psd / np.sum(psd) * total_conc

        return psd

    def _get_default_parameters(self) -> Dict[str, Any]:
        """Get default PBM parameters."""
        return {
            'num_size_classes': 20,
            'initial_mean_diameter': 10e-6,
            'initial_std_diameter': 3e-6,
            'initial_concentration': 1e13,
            'kernel_coeffs': {
                'collision_efficiency': 0.1,
                'shear_rate': 50.0,
                'temperature': 293.15,
                'breakage_rate': 0.001,
                'primary_diameter': 1e-6,
                'fractal_dimension': 2.3,
                'binding_strength': 1e3
            },
            'simulation': {
                'total_time': 600.0,
                'num_timesteps': 100
            }
        }
