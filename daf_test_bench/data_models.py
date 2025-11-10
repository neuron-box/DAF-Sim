"""
Data Models for DAF Test Bench

Defines standardized data structures for test configurations and benchmark results
as specified in Section 4.0 of DAF-TB-HLD-v1.0.

These models ensure consistent data interchange between the Test Harness,
engine wrappers, and metrics collection system (FR-4.1, FR-5).
"""

from dataclasses import dataclass, field, asdict
from typing import Dict, Any, Optional, List
from datetime import datetime
import json


@dataclass
class HydrodynamicStimuli:
    """Hydrodynamic input parameters."""
    inlet_flow_rate: float  # [m³/s]
    recycle_ratio: float  # [-]
    operating_pressure: float  # [Pa]
    tank_geometry: Optional[Dict[str, float]] = None
    additional_params: Dict[str, Any] = field(default_factory=dict)


@dataclass
class PBMStimuli:
    """Population Balance Model input parameters."""
    influent_psd: List[float]  # Particle size distribution [#/m³]
    kernel_coeffs: Dict[str, float]  # Aggregation/breakage coefficients
    size_classes: Optional[List[float]] = None  # Size class boundaries [m]
    additional_params: Dict[str, Any] = field(default_factory=dict)


@dataclass
class SolverStimuli:
    """Solver configuration parameters."""
    timestep: float  # [s]
    convergence_limit: float  # [-]
    max_iterations: Optional[int] = 10000
    relaxation_factors: Optional[Dict[str, float]] = None
    additional_params: Dict[str, Any] = field(default_factory=dict)


@dataclass
class TestStimuli:
    """Collection of all standardized stimuli for a test case."""
    hydrodynamics: Optional[HydrodynamicStimuli] = None
    pbm: Optional[PBMStimuli] = None
    solver: Optional[SolverStimuli] = None
    additional: Dict[str, Any] = field(default_factory=dict)


@dataclass
class TestConfiguration:
    """
    Standardized test configuration (config.json schema).

    This data model defines the complete specification for a test case,
    including all input parameters (stimuli), mesh files, and ground truth
    validation data paths (FR-4.1).
    """
    test_name: str
    description: str
    mesh_file: Optional[str] = None
    stimuli: Optional[TestStimuli] = None
    ground_truth_files: Dict[str, str] = field(default_factory=dict)
    engine_specific_settings: Dict[str, Dict[str, Any]] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'TestConfiguration':
        """
        Create TestConfiguration from dictionary (parsed JSON).

        Args:
            data: Dictionary containing test configuration

        Returns:
            TestConfiguration instance
        """
        # Parse stimuli
        stimuli_data = data.get('stimuli', {})
        stimuli = None
        if stimuli_data:
            hydro = None
            if 'hydrodynamics' in stimuli_data:
                hydro = HydrodynamicStimuli(**stimuli_data['hydrodynamics'])

            pbm = None
            if 'pbm' in stimuli_data:
                pbm = PBMStimuli(**stimuli_data['pbm'])

            solver = None
            if 'solver' in stimuli_data:
                solver = SolverStimuli(**stimuli_data['solver'])

            stimuli = TestStimuli(
                hydrodynamics=hydro,
                pbm=pbm,
                solver=solver,
                additional=stimuli_data.get('additional', {})
            )

        return cls(
            test_name=data['test_name'],
            description=data['description'],
            mesh_file=data.get('mesh_file'),
            stimuli=stimuli,
            ground_truth_files=data.get('ground_truth_files', {}),
            engine_specific_settings=data.get('engine_specific_settings', {})
        )

    @classmethod
    def from_json_file(cls, filepath: str) -> 'TestConfiguration':
        """
        Load TestConfiguration from JSON file.

        Args:
            filepath: Path to config.json file

        Returns:
            TestConfiguration instance
        """
        with open(filepath, 'r') as f:
            data = json.load(f)
        return cls.from_dict(data)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return asdict(self)


@dataclass
class ScientificMetrics:
    """
    Scientific performance metrics (FR-5.1).

    All metrics related to physical accuracy and model validation.
    """
    particle_removal_eff: Optional[float] = None  # [0-1]
    velocity_error_l2: Optional[float] = None  # L2 norm
    tke_error_l2: Optional[float] = None  # L2 norm
    outlet_psd_error_wasserstein: Optional[float] = None  # Wasserstein distance
    mean_floc_size_d50: Optional[float] = None  # [m]
    mean_residence_time: Optional[float] = None  # [s]
    mass_conservation_error_pct: Optional[float] = None  # [%]
    additional_metrics: Dict[str, float] = field(default_factory=dict)


@dataclass
class ComputationalMetrics:
    """
    Computational performance metrics (FR-5.2).

    All metrics related to computational efficiency and resource usage.
    """
    wall_clock_sec: Optional[float] = None  # [s]
    cpu_hours: Optional[float] = None  # [hours]
    peak_ram_gb: Optional[float] = None  # [GB]
    converged: Optional[bool] = None
    num_iterations: Optional[int] = None
    scalability_metrics: Dict[str, float] = field(default_factory=dict)


@dataclass
class BenchmarkResult:
    """
    Complete benchmark result for a single test run.

    This data model defines the standardized output format (results.json schema)
    that each engine must produce (Section 4.0, DAF-TB-HLD-v1.0).
    """
    engine_name: str
    test_name: str
    timestamp: str
    run_status: str  # "Success" or "Failed"
    scientific_metrics: ScientificMetrics
    computational_metrics: ComputationalMetrics
    run_log: str = ""

    @classmethod
    def create(
        cls,
        engine_name: str,
        test_name: str,
        run_status: str = "Success",
        sci_metrics: Optional[ScientificMetrics] = None,
        comp_metrics: Optional[ComputationalMetrics] = None,
        run_log: str = ""
    ) -> 'BenchmarkResult':
        """
        Create a BenchmarkResult with automatic timestamp.

        Args:
            engine_name: Name of the engine
            test_name: Name of the test case
            run_status: "Success" or "Failed"
            sci_metrics: Scientific metrics (defaults to empty)
            comp_metrics: Computational metrics (defaults to empty)
            run_log: Execution log

        Returns:
            BenchmarkResult instance
        """
        return cls(
            engine_name=engine_name,
            test_name=test_name,
            timestamp=datetime.utcnow().isoformat() + "Z",
            run_status=run_status,
            scientific_metrics=sci_metrics or ScientificMetrics(),
            computational_metrics=comp_metrics or ComputationalMetrics(),
            run_log=run_log
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return asdict(self)

    def to_json(self, indent: int = 2) -> str:
        """Convert to JSON string."""
        return json.dumps(self.to_dict(), indent=indent)

    def save_to_file(self, filepath: str) -> None:
        """
        Save result to JSON file.

        Args:
            filepath: Path to output results.json file
        """
        with open(filepath, 'w') as f:
            f.write(self.to_json())

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'BenchmarkResult':
        """Create BenchmarkResult from dictionary."""
        return cls(
            engine_name=data['engine_name'],
            test_name=data['test_name'],
            timestamp=data['timestamp'],
            run_status=data['run_status'],
            scientific_metrics=ScientificMetrics(**data['scientific_metrics']),
            computational_metrics=ComputationalMetrics(**data['computational_metrics']),
            run_log=data.get('run_log', '')
        )

    @classmethod
    def from_json_file(cls, filepath: str) -> 'BenchmarkResult':
        """Load BenchmarkResult from JSON file."""
        with open(filepath, 'r') as f:
            data = json.load(f)
        return cls.from_dict(data)
