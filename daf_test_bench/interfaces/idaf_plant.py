"""
IDAFPlant Interface

Abstract Base Class defining the standard interface for all DAF simulation engines.
Implements the pluggable architecture per FR-1.1 (DAF-TB-SRR-v1.0).

This interface enables the Test Harness to interact with heterogeneous DAF engines
(Python, C++, OpenFOAM, etc.) through a unified API, following the Abstract Base
Class pattern [39, 40, 41, 42, 43].
"""

from abc import ABC, abstractmethod
from typing import Dict, Any, Optional
import numpy as np


class IDAFPlant(ABC):
    """
    Abstract Base Class defining the standard interface for all DAF simulation engines.

    All DAF Plant engines must implement this interface to be compatible with the
    Test Harness. This follows the Component-Based Software Engineering (CBSE)
    design pattern to promote separation of concerns (NFR-1).

    Workflow:
        1. setup(config) - Configure simulation with standardized stimuli
        2. initialize() - Initialize simulation state, load mesh, etc.
        3. run() - Execute simulation to completion
        4. get_metrics() - Retrieve performance metrics
        5. get_field_data(var) - Retrieve field data for validation
        6. finalize() - Clean up resources
    """

    def __init__(self, engine_name: str):
        """
        Initialize the DAF Plant engine.

        Args:
            engine_name: Unique identifier for this engine (e.g., "Pillar3_Physics")
        """
        self.engine_name = engine_name
        self._is_initialized = False
        self._is_setup = False
        self._run_completed = False
        self._config: Optional[Dict[str, Any]] = None

    @abstractmethod
    def setup(self, configuration: Dict[str, Any]) -> bool:
        """
        Configure the simulation with standardized stimuli.

        This method receives the TestConfiguration dictionary (parsed from config.json)
        and translates it into engine-specific configuration files or data structures.

        Args:
            configuration: Standardized test configuration dictionary containing:
                - test_name: Name of the test case
                - description: Test description
                - mesh_file: Path to mesh file (if applicable)
                - stimuli: Dictionary of input parameters
                    - hydrodynamics: Flow parameters
                    - pbm: Population balance parameters
                    - solver: Solver settings
                - ground_truth_files: Paths to experimental validation data
                - engine_specific_settings: Engine-specific parameters

        Returns:
            bool: True if setup successful, False otherwise

        Raises:
            ValueError: If configuration is invalid
        """
        pass

    @abstractmethod
    def initialize(self) -> bool:
        """
        Initialize the simulation state.

        This method performs all pre-processing steps required before running:
        - Loading mesh/grid
        - Initializing field variables
        - Setting initial conditions
        - Allocating memory

        Must be called after setup() and before run().

        Returns:
            bool: True if initialization successful, False otherwise
        """
        pass

    @abstractmethod
    def run(self) -> bool:
        """
        Run the entire simulation to completion.

        This method executes the main simulation loop until:
        - Steady-state is reached (for steady problems)
        - Final time is reached (for transient problems)
        - Convergence criteria are met
        - Maximum iterations exceeded

        Must be called after initialize().

        Returns:
            bool: True if simulation completed successfully, False if failed
        """
        pass

    @abstractmethod
    def get_metrics(self) -> Dict[str, Any]:
        """
        Return all performance metrics in standardized format.

        This method extracts and returns both scientific and computational metrics
        as defined in FR-5.1 and FR-5.2 (DAF-TB-SRR-v1.0).

        Must be called after run().

        Returns:
            dict: BenchmarkResult dictionary containing:
                - engine_name: Name of this engine
                - test_name: Name of the test case
                - timestamp: ISO 8601 timestamp
                - run_status: "Success" or "Failed"
                - scientific_metrics: Dict of scientific KPIs
                    - particle_removal_eff: Removal efficiency [0-1]
                    - velocity_error_l2: L2 norm of velocity error
                    - tke_error_l2: L2 norm of TKE error
                    - outlet_psd_error_wasserstein: Wasserstein distance for PSD
                    - mean_floc_size_d50: Mean floc diameter [m]
                    - mean_residence_time: Mean residence time [s]
                    - mass_conservation_error_pct: Mass balance error [%]
                - computational_metrics: Dict of computational KPIs
                    - wall_clock_sec: Wall-clock time [s]
                    - cpu_hours: CPU time [hours]
                    - peak_ram_gb: Peak memory usage [GB]
                    - converged: Whether simulation converged (bool)
                - run_log: Full stdout/stderr log
        """
        pass

    @abstractmethod
    def get_field_data(self, variable_name: str) -> Optional[np.ndarray]:
        """
        Return full field data for validation.

        This method extracts complete field variables (velocity, pressure, PSD, etc.)
        for comparison against experimental ground truth data during V&V.

        Args:
            variable_name: Name of field variable to extract (e.g., "velocity",
                          "pressure", "tke", "psd", "concentration")

        Returns:
            np.ndarray: Field data array, or None if variable not available
                       Shape depends on problem dimension and variable type
        """
        pass

    @abstractmethod
    def finalize(self) -> None:
        """
        Clean up all resources.

        This method performs all cleanup operations:
        - Closing files
        - Freeing memory
        - Terminating subprocesses
        - Deleting temporary files

        Should be called after all metrics and field data have been extracted.
        """
        pass

    # Helper properties for state tracking
    @property
    def is_setup(self) -> bool:
        """Check if engine has been configured."""
        return self._is_setup

    @property
    def is_initialized(self) -> bool:
        """Check if engine has been initialized."""
        return self._is_initialized

    @property
    def run_completed(self) -> bool:
        """Check if simulation run has completed."""
        return self._run_completed

    def validate_workflow(self, step: str) -> None:
        """
        Validate that workflow steps are called in correct order.

        Args:
            step: Current workflow step being executed

        Raises:
            RuntimeError: If workflow order is violated
        """
        if step == "initialize" and not self._is_setup:
            raise RuntimeError("Cannot initialize(): Must call setup() first")
        if step == "run" and not self._is_initialized:
            raise RuntimeError("Cannot run(): Must call initialize() first")
        if step in ["get_metrics", "get_field_data"] and not self._run_completed:
            raise RuntimeError(f"Cannot {step}(): Must call run() first")
