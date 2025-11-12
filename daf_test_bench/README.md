# DAF Test Bench Framework

## Overview

The **DAF Test Bench** is a unified testing and benchmarking framework for Dissolved Air Flotation (DAF) simulation engines. It implements a pluggable architecture that enables multiple heterogeneous simulation engines to be tested, validated, and compared using standardized test cases and metrics.

This framework implements the specifications defined in:
- **DAF-TB-SRR-v1.0**: Software Requirements Review
- **DAF-TB-HLD-v1.0**: High-Level Design

## Architecture

The framework follows a **Component-Based Software Engineering (CBSE)** design with three key patterns:

1. **Abstract Base Class (ABC) Pattern**: The `IDAFPlant` interface defines the standard contract that all DAF engines must implement.

2. **Facade Pattern**: Engine-specific wrappers encapsulate implementation details and provide a uniform interface to the Test Harness.

3. **Strategy Pattern**: Different engines can be plugged in and executed interchangeably through the common interface.

```
┌─────────────────────┐
│   Test Harness      │
│  (Orchestrator)     │
└──────────┬──────────┘
           │
           │ uses
           ▼
┌──────────────────────┐
│  <<Interface>>       │
│    IDAFPlant         │
├──────────────────────┤
│ + setup()            │
│ + initialize()       │
│ + run()              │
│ + get_metrics()      │
│ + get_field_data()   │
│ + finalize()         │
└──────────┬───────────┘
           │
           │ implemented by
           ▼
┌────────────────────────────────────────────┐
│  Engine Wrappers (Facades)                 │
├────────────────┬───────────────┬───────────┤
│ Pillar3Physics │ Floc_Kinetics │ ...More..│
│   Wrapper      │  PBM_Wrapper  │           │
└────────────────┴───────────────┴───────────┘
```

## Features

✅ **Pluggable Architecture**: Add new engines without modifying core framework
✅ **Standardized I/O**: Common `config.json` input and `results.json` output format
✅ **Comprehensive Metrics**: Scientific and computational performance tracking
✅ **Automated Benchmarking**: Run test suites across multiple engines
✅ **Export Capabilities**: JSON, CSV, and text report formats
✅ **Workflow Validation**: Enforces correct execution order
✅ **Test Suite**: Comprehensive unit and integration tests

## Installation

### Requirements

- Python 3.8+
- NumPy >= 1.22.0
- SciPy >= 1.10.0
- pandas >= 1.3.0
- psutil >= 5.9.0

### Install the Framework

```bash
cd DAF-Sim
pip install -e .
```

This installs the `daf-test-bench` package and all its dependencies.

### Install DAF Engines

The framework requires individual DAF engine components to be installed:

```bash
# Install Pillar3 Physics Model
cd pillar3_physics_model
pip install -e .

# Install Floc Kinetics PBM
cd ../floc_kinetics_pbm
pip install -e .
```

## Quick Start

### Running a Single Benchmark

```bash
# Run a test case on all available engines
python daf_benchmark.py --config daf_test_bench/test_cases/T3_PBM/config.json

# Run on a specific engine
python daf_benchmark.py \
  --config daf_test_bench/test_cases/T2_Multiphase/config.json \
  --engines Pillar3_Physics_Model
```

### Running a Test Suite

```bash
# Run all test cases
python daf_benchmark.py \
  --suite daf_test_bench/test_cases/ \
  --output results/

# Run only T3 tests
python daf_benchmark.py \
  --suite daf_test_bench/test_cases/ \
  --filter T3_ \
  --output results/
```

### List Available Engines

```bash
python daf_benchmark.py --list-engines
```

## Usage Examples

### Example 1: Benchmarking Pillar3 Physics Model

```python
from daf_test_bench import TestHarness
from daf_test_bench.engines import Pillar3PhysicsWrapper

# Create test harness
harness = TestHarness()

# Register engine
harness.register_engine(Pillar3PhysicsWrapper)

# Run benchmark
collector = harness.run_benchmark(
    config_file="daf_test_bench/test_cases/T2_Multiphase/config.json",
    verbose=True
)

# Save results
harness.save_results("results/")
```

### Example 2: Programmatic Engine Execution

```python
from daf_test_bench.engines import Pillar3PhysicsWrapper

# Create engine instance
engine = Pillar3PhysicsWrapper()

# Configure
config = {
    'test_name': 'Custom_Test',
    'description': 'Testing collision efficiency',
    'stimuli': {
        'pillar3': {
            'particle': {
                'diameter': 15e-6,
                'zeta_potential': -0.010,
                'density': 1200.0,
                'hamaker_constant': 8e-21
            },
            'bubble': {
                'diameter': 40e-6,
                'zeta_potential': -0.020
            },
            'fluid': {
                'use_water_properties': True,
                'ionic_strength': 0.01,
                'ph': 7.0
            }
        }
    }
}

# Execute workflow
engine.setup(config)
engine.initialize()
engine.run()

# Get results
metrics = engine.get_metrics()
print(f"Collision efficiency: {metrics['scientific_metrics']['additional_metrics']['collision_efficiency_alpha_bp']}")

engine.finalize()
```

## Configuration Format

### Input: config.json

```json
{
  "test_name": "T3_PBM_Jar_Test",
  "description": "Population Balance Model test case",
  "mesh_file": null,
  "stimuli": {
    "pbm": {
      "kernel_coeffs": {
        "collision_efficiency": 0.15,
        "shear_rate": 75.0
      },
      "simulation": {
        "total_time": 900.0,
        "num_timesteps": 150
      }
    },
    "solver": {
      "timestep": 6.0,
      "convergence_limit": 1e-6
    }
  },
  "ground_truth_files": {
    "fbrm": "path/to/experimental_data.csv"
  },
  "engine_specific_settings": {
    "Engine_Name": {
      "custom_parameter": "value"
    }
  }
}
```

### Output: results.json

```json
{
  "engine_name": "Pillar3_Physics_Model",
  "test_name": "T2_Multiphase",
  "timestamp": "2025-11-10T21:30:00Z",
  "run_status": "Success",
  "scientific_metrics": {
    "particle_removal_eff": 0.95,
    "velocity_error_l2": 0.045,
    "mean_floc_size_d50": 45.5e-6,
    "mass_conservation_error_pct": 0.001
  },
  "computational_metrics": {
    "wall_clock_sec": 120.5,
    "cpu_hours": 0.0335,
    "peak_ram_gb": 0.5,
    "converged": true
  },
  "run_log": "Full execution log..."
}
```

## Test Cases

The framework includes several standardized test cases:

- **T1_Hydro**: Single-phase hydrodynamics (steady-state)
- **T2_Multiphase**: Two-phase flow with bubble-particle interactions
- **T3_PBM**: 0D Population Balance Model (jar test)
- **T4_Transport**: Coupled transport with passive scalars
- **T5_FullSystem**: Fully coupled DAF simulation

## Implementing a New Engine

To integrate a new DAF engine:

### 1. Create a Wrapper Class

```python
from daf_test_bench.interfaces.idaf_plant import IDAFPlant
from daf_test_bench.data_models import BenchmarkResult, ScientificMetrics, ComputationalMetrics
import numpy as np

class MyEngineWrapper(IDAFPlant):
    def __init__(self):
        super().__init__(engine_name="My_Engine")
        # Initialize your engine-specific state

    def setup(self, configuration):
        # Parse config and set up your engine
        self._is_setup = True
        return True

    def initialize(self):
        self.validate_workflow("initialize")
        # Initialize simulation state
        self._is_initialized = True
        return True

    def run(self):
        self.validate_workflow("run")
        # Run your simulation
        self._run_completed = True
        return True

    def get_metrics(self):
        self.validate_workflow("get_metrics")
        # Extract metrics in standardized format
        sci_metrics = ScientificMetrics(particle_removal_eff=0.95)
        comp_metrics = ComputationalMetrics(wall_clock_sec=100.0, converged=True)

        result = BenchmarkResult.create(
            engine_name=self.engine_name,
            test_name=self._config['test_name'],
            run_status="Success",
            sci_metrics=sci_metrics,
            comp_metrics=comp_metrics
        )
        return result.to_dict()

    def get_field_data(self, variable_name):
        # Return field data arrays
        return None

    def finalize(self):
        # Clean up resources
        pass
```

### 2. Register the Engine

```python
from daf_test_bench import TestHarness
from my_engine_wrapper import MyEngineWrapper

harness = TestHarness()
harness.register_engine(MyEngineWrapper)
```

## API Reference

### IDAFPlant Interface

**Core Methods:**
- `setup(configuration: dict) -> bool`: Configure simulation
- `initialize() -> bool`: Initialize state
- `run() -> bool`: Execute simulation
- `get_metrics() -> dict`: Extract performance metrics
- `get_field_data(variable_name: str) -> np.ndarray`: Get field data
- `finalize() -> None`: Clean up resources

**Properties:**
- `engine_name: str`: Unique engine identifier
- `is_setup: bool`: Setup completion status
- `is_initialized: bool`: Initialization status
- `run_completed: bool`: Run completion status

### TestHarness

**Methods:**
- `register_engine(engine_class)`: Register a DAF engine
- `list_engines() -> List[str]`: Get registered engine names
- `run_benchmark(config_file, engines, verbose)`: Run single benchmark
- `run_test_suite(test_cases_dir, engines, test_filter)`: Run test suite
- `save_results(output_dir)`: Export all results

### MetricsCollector

**Methods:**
- `add_result(result)`: Add a benchmark result
- `get_results_for_test(test_name)`: Filter by test
- `get_results_for_engine(engine_name)`: Filter by engine
- `export_to_json(filepath)`: Export to JSON
- `export_to_csv(filepath)`: Export to CSV
- `generate_report(output_dir)`: Generate complete report

## Testing

Run the comprehensive test suite:

```bash
python daf_test_bench/tests/test_framework.py
```

## Project Structure

```
daf_test_bench/
├── __init__.py                    # Package initialization
├── interfaces/
│   ├── __init__.py
│   └── idaf_plant.py             # IDAFPlant ABC
├── engines/
│   ├── __init__.py
│   ├── pillar3_wrapper.py        # Pillar3 Physics wrapper
│   └── floc_pbm_wrapper.py       # Floc PBM wrapper
├── metrics/
│   ├── __init__.py
│   └── metrics_collector.py      # Metrics collection
├── data_models.py                # Data structures
├── test_harness.py               # Main orchestrator
├── test_cases/                   # Test configurations
│   ├── T1_Hydro/
│   ├── T2_Multiphase/
│   ├── T3_PBM/
│   ├── T4_Transport/
│   └── T5_FullSystem/
├── tests/
│   ├── __init__.py
│   └── test_framework.py         # Test suite
├── requirements.txt              # Dependencies
└── README.md                     # This file
```

## Contributing

To add a new DAF engine:

1. Implement the `IDAFPlant` interface
2. Create a wrapper in `engines/`
3. Add engine registration to `daf_benchmark.py`
4. Create test cases in `test_cases/`
5. Add unit tests in `tests/`

## References

- **DAF-TB-SRR-v1.0**: Software Requirements Review
- **DAF-TB-HLD-v1.0**: High-Level Design Document
- Component-Based Software Engineering (CBSE) principles
- Abstract Base Class pattern (Python ABC module)
- Facade design pattern

## License

See LICENSE file for details.

## Authors

Developed for the DAF-Sim multi-pillar simulator project.

## Version

Framework Version: 1.0.0
