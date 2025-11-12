# DAF Test Bench - Quick Start Guide

## Installation

### 1. Install the DAF Test Bench Package

```bash
cd DAF-Sim
pip install -e .
```

This will automatically install all required dependencies (numpy, scipy, pandas, psutil).

### 2. Install DAF Engine Components

```bash
# Install Pillar3 Physics Model
cd pillar3_physics_model
pip install -e .
cd ..

# Install Floc Kinetics PBM
cd floc_kinetics_pbm
pip install -e .
cd ..
```

## Running Your First Benchmark

### List Available Engines

```bash
python daf_benchmark.py --list-engines
```

Expected output:
```
Available DAF Plant Engines:
  - Pillar3_Physics_Model
  - Floc_Kinetics_PBM
```

### Run a Single Test Case

```bash
# Run T2 (Multiphase) test on Pillar3 Physics Model
python daf_benchmark.py \
  --config daf_test_bench/test_cases/T2_Multiphase/config.json \
  --engines Pillar3_Physics_Model \
  --output results/T2_test
```

### Run a Test Suite

```bash
# Run all available test cases on all engines
python daf_benchmark.py \
  --suite daf_test_bench/test_cases/ \
  --output results/full_suite
```

### View Results

After running, check the output directory:

```bash
ls -la results/
cat results/report.txt           # Human-readable summary
cat results/full_results.json    # Complete JSON data
cat results/results_table.csv    # Tabular format for Excel/plotting
```

## Example Output

### Console Output

```
================================================================================
DAF Test Bench - Starting Benchmark
================================================================================
Loading configuration: daf_test_bench/test_cases/T2_Multiphase/config.json
Test: T2_Multiphase_Standard_k-Epsilon
Description: Standard benchmark for two-phase hydrodynamics...

--------------------------------------------------------------------------------
Running: Pillar3_Physics_Model
--------------------------------------------------------------------------------
  1. Setting up...
  2. Initializing...
  3. Running simulation...
  4. Collecting metrics...
  5. Finalizing...
Status: Success
Wall Clock: 0.003 s

================================================================================
Benchmark Complete
================================================================================
Total runs: 1
Successful: 1/1
```

### Results JSON

```json
{
  "engine_name": "Pillar3_Physics_Model",
  "test_name": "T2_Multiphase_Standard_k-Epsilon",
  "timestamp": "2025-11-10T21:45:00Z",
  "run_status": "Success",
  "scientific_metrics": {
    "additional_metrics": {
      "collision_efficiency_alpha_bp": 0.2345,
      "eta_collision": 0.4567,
      "alpha_attachment": 0.5134
    }
  },
  "computational_metrics": {
    "wall_clock_sec": 0.003,
    "cpu_hours": 8.33e-07,
    "peak_ram_gb": 0.1,
    "converged": true
  }
}
```

## Python API Usage

### Simple Example

```python
from daf_test_bench import TestHarness
from daf_test_bench.engines import Pillar3PhysicsWrapper

# Create and configure harness
harness = TestHarness()
harness.register_engine(Pillar3PhysicsWrapper)

# Run benchmark
results = harness.run_benchmark(
    config_file="daf_test_bench/test_cases/T2_Multiphase/config.json",
    verbose=True
)

# Save results
harness.save_results("my_results/")
```

### Direct Engine Interaction

```python
from daf_test_bench.engines import Pillar3PhysicsWrapper

engine = Pillar3PhysicsWrapper()

# Minimal configuration
config = {
    'test_name': 'Quick_Test',
    'description': 'Quick collision efficiency test',
    'stimuli': {
        'pillar3': {
            'particle': {'diameter': 15e-6, 'zeta_potential': -0.01,
                        'density': 1200.0, 'hamaker_constant': 8e-21},
            'bubble': {'diameter': 40e-6, 'zeta_potential': -0.02},
            'fluid': {'use_water_properties': True, 'ionic_strength': 0.01, 'ph': 7.0}
        }
    }
}

# Execute
engine.setup(config)
engine.initialize()
engine.run()

# Get results
metrics = engine.get_metrics()
alpha_bp = metrics['scientific_metrics']['additional_metrics']['collision_efficiency_alpha_bp']
print(f"Collision Efficiency (α_bp): {alpha_bp:.4f}")

# Get field data
efficiency_data = engine.get_field_data('alpha_bp')
print(f"Raw value: {efficiency_data[0]:.6f}")

engine.finalize()
```

## Creating Your Own Test Case

### 1. Create Test Directory

```bash
mkdir -p daf_test_bench/test_cases/T_Custom
```

### 2. Create config.json

```json
{
  "test_name": "T_Custom_MyTest",
  "description": "Custom test case description",
  "stimuli": {
    "pillar3": {
      "particle": {
        "diameter": 20e-6,
        "zeta_potential": -0.015,
        "density": 1300.0,
        "hamaker_constant": 9e-21
      },
      "bubble": {
        "diameter": 50e-6,
        "zeta_potential": -0.025
      },
      "fluid": {
        "use_water_properties": true,
        "ionic_strength": 0.012,
        "ph": 7.5
      }
    }
  }
}
```

### 3. Run Your Test

```bash
python daf_benchmark.py \
  --config daf_test_bench/test_cases/T_Custom/config.json \
  --output results/custom_test
```

## Troubleshooting

### Issue: "No module named 'numpy'"

**Solution**: Install dependencies
```bash
pip install numpy scipy pandas
```

### Issue: "Cannot import pillar3_physics_model"

**Solution**: Install engine packages
```bash
cd pillar3_physics_model && pip install -e .
cd ../floc_kinetics_pbm && pip install -e .
```

### Issue: Empty results or "Failed" status

**Solution**: Check the run log in results.json:
```python
import json
with open('results/full_results.json') as f:
    data = json.load(f)
    print(data['results'][0]['run_log'])
```

## Next Steps

1. **Explore the examples**: See `daf_test_bench/test_cases/` for more configurations
2. **Read the API docs**: See `daf_test_bench/README.md`
3. **Create custom engines**: Implement the `IDAFPlant` interface for your own simulators
4. **Run validation tests**: Use the test suite in `daf_test_bench/tests/`

## Support

For issues or questions:
- Check `daf_test_bench/README.md` for detailed documentation
- Review test cases in `daf_test_bench/test_cases/`
- See example usage in `daf_test_bench/tests/test_framework.py`
