"""
Pillar3 Physics Model Engine Wrapper

Facade/Adapter for the Pillar3 Physics Model (collision efficiency calculation).
Implements the IDAFPlant interface for integration with the Test Harness.

This wrapper translates standardized test configurations into Pillar3-specific
API calls and extracts results in the standardized BenchmarkResult format.
"""

import time
import traceback
import logging
from typing import Dict, Any, Optional
import numpy as np

from pillar3_physics_model import (
    Particle, Bubble, FluidProperties,
    calculate_attachment_efficiency
)

from ..interfaces.idaf_plant import IDAFPlant
from ..data_models import BenchmarkResult, ScientificMetrics, ComputationalMetrics


class Pillar3PhysicsWrapper(IDAFPlant):
    """
    Wrapper for Pillar3 Physics Model engine.

    The Pillar3 Physics Model calculates collision efficiency factors (α_bp)
    for bubble-particle interactions using trajectory analysis and DLVO theory.
    This is typically used as a sub-model within larger DAF simulations.
    """

    def __init__(self):
        """Initialize Pillar3 Physics Model wrapper."""
        super().__init__(engine_name="Pillar3_Physics_Model")
        self._particle: Optional[Particle] = None
        self._bubble: Optional[Bubble] = None
        self._fluid: Optional[FluidProperties] = None
        self._results: Optional[Dict[str, Any]] = None
        self._start_time: float = 0.0
        self._end_time: float = 0.0

        # Set up logging
        self.logger = logging.getLogger(f"{__name__}.{self.engine_name}")

    def setup(self, configuration: Dict[str, Any]) -> bool:
        """
        Configure the Pillar3 Physics Model simulation.

        Expected configuration structure:
            stimuli:
                pillar3:
                    particle:
                        diameter: [m]
                        zeta_potential: [V]
                        density: [kg/m³]
                        hamaker_constant: [J]
                    bubble:
                        diameter: [m]
                        zeta_potential: [V]
                    fluid:
                        temperature: [K]
                        ionic_strength: [M]
                        ph: [-]
                        viscosity: [Pa·s]
                        density: [kg/m³]

        Args:
            configuration: Test configuration dictionary

        Returns:
            bool: True if setup successful
        """
        try:
            self._config = configuration
            self.logger.info(f"Setting up {self.engine_name}")
            self.logger.info(f"Test: {configuration.get('test_name', 'Unknown')}")

            # Extract engine-specific settings
            engine_settings = configuration.get('engine_specific_settings', {}).get(
                self.engine_name, {}
            )

            # Extract Pillar3 stimuli
            stimuli = configuration.get('stimuli', {})
            pillar3_params = stimuli.get('pillar3', {})

            if not pillar3_params:
                # Use default DAF parameters if not specified
                pillar3_params = self._get_default_parameters()
                self.logger.info("Using default DAF parameters")

            # Create particle
            particle_params = pillar3_params.get('particle', {})
            self._particle = Particle(
                diameter=particle_params.get('diameter', 15e-6),
                zeta_potential=particle_params.get('zeta_potential', -0.010),
                density=particle_params.get('density', 1200.0),
                hamaker_constant=particle_params.get('hamaker_constant', 8e-21)
            )

            # Create bubble
            bubble_params = pillar3_params.get('bubble', {})
            self._bubble = Bubble(
                diameter=bubble_params.get('diameter', 40e-6),
                zeta_potential=bubble_params.get('zeta_potential', -0.020)
            )

            # Create fluid properties
            fluid_params = pillar3_params.get('fluid', {})
            if fluid_params.get('use_water_properties', True):
                self._fluid = FluidProperties.water_at_20C(
                    ionic_strength=fluid_params.get('ionic_strength', 0.01),
                    ph=fluid_params.get('ph', 7.0)
                )
            else:
                # Custom fluid properties
                self._fluid = FluidProperties(
                    temperature=fluid_params.get('temperature', 293.15),
                    viscosity=fluid_params.get('viscosity', 1.002e-3),
                    density=fluid_params.get('density', 998.2),
                    ionic_strength=fluid_params.get('ionic_strength', 0.01),
                    relative_permittivity=fluid_params.get('relative_permittivity', 78.5)
                )

            self.logger.info(f"Particle diameter: {self._particle.diameter*1e6:.1f} μm")
            self.logger.info(f"Bubble diameter: {self._bubble.diameter*1e6:.1f} μm")
            self.logger.info(f"Ionic strength: {self._fluid.ionic_strength:.4f} M")

            self._is_setup = True
            return True

        except Exception as e:
            self.logger.error(f"ERROR in setup: {str(e)}")
            self.logger.debug(traceback.format_exc())
            return False

    def initialize(self) -> bool:
        """
        Initialize the simulation.

        For Pillar3, this is a simple physics model with no mesh or field
        initialization required.

        Returns:
            bool: True if initialization successful
        """
        try:
            self.validate_workflow("initialize")
            self.logger.info("Initializing Pillar3 Physics Model")

            # Verify all components are configured
            if self._particle is None or self._bubble is None or self._fluid is None:
                raise RuntimeError("Not all components configured")

            self._is_initialized = True
            self.logger.info("Initialization complete")
            return True

        except Exception as e:
            self.logger.error(f"ERROR in initialize: {str(e)}")
            self.logger.debug(traceback.format_exc())
            return False

    def run(self) -> bool:
        """
        Run the collision efficiency calculation.

        Returns:
            bool: True if calculation successful
        """
        try:
            self.validate_workflow("run")
            self.logger.info("Running collision efficiency calculation")

            self._start_time = time.time()

            # Calculate collision efficiency using trajectory analysis + DLVO
            self._results = calculate_attachment_efficiency(
                self._particle,
                self._bubble,
                self._fluid
            )

            self._end_time = time.time()
            elapsed = self._end_time - self._start_time

            self.logger.info(f"Calculation complete in {elapsed:.6f} s")
            self.logger.info(f"α_bp = {self._results['alpha_bp']:.6f}")
            self.logger.info(f"η_collision = {self._results['eta_collision']:.6f}")
            self.logger.info(f"α_attachment = {self._results['alpha_attachment']:.6f}")

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
            if self._results:
                # Store collision efficiency as a custom metric
                sci_metrics.additional_metrics = {
                    'collision_efficiency_alpha_bp': self._results['alpha_bp'],
                    'eta_collision': self._results['eta_collision'],
                    'alpha_attachment': self._results['alpha_attachment'],
                    'critical_offset_m': self._results['critical_offset'],
                    'debye_length_m': self._results['debye_length'],
                    'energy_barrier_kT': self._results['energy_barrier_kT']
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
                num_iterations=1  # Single calculation
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

        Pillar3 is a point model (not spatially resolved), so field data
        is limited to scalar results.

        Args:
            variable_name: Name of variable (e.g., "alpha_bp", "energy_barrier")

        Returns:
            Scalar value as 0D array, or None if not available
        """
        if not self._results:
            return None

        mapping = {
            'alpha_bp': 'alpha_bp',
            'collision_efficiency': 'alpha_bp',
            'eta_collision': 'eta_collision',
            'alpha_attachment': 'alpha_attachment',
            'critical_offset': 'critical_offset',
            'debye_length': 'debye_length',
            'energy_barrier': 'energy_barrier_kT'
        }

        key = mapping.get(variable_name)
        if key and key in self._results:
            return np.array([self._results[key]])

        return None

    def finalize(self) -> None:
        """Clean up resources."""
        self.logger.info("Finalizing Pillar3 Physics Model")
        # No cleanup needed for this simple model

    def _get_log_contents(self) -> str:
        """
        Get logging contents.

        In the future, this could return contents from a logging handler.
        For now, returns a simple summary.
        """
        return f"Pillar3 Physics Model log for {self.engine_name}"

    def _get_default_parameters(self) -> Dict[str, Any]:
        """Get default DAF operating parameters."""
        return {
            'particle': {
                'diameter': 15e-6,  # 15 μm (flocculated)
                'zeta_potential': -0.010,  # -10 mV (after coagulation)
                'density': 1200.0,  # kg/m³
                'hamaker_constant': 8e-21  # J (hydrophobic)
            },
            'bubble': {
                'diameter': 40e-6,  # 40 μm (typical DAF)
                'zeta_potential': -0.020  # -20 mV
            },
            'fluid': {
                'use_water_properties': True,
                'ionic_strength': 0.01,  # 0.01 M (after coagulant)
                'ph': 7.0
            }
        }
