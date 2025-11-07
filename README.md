# DAF-Sim
Hub for Dissolved Air Flotation (DAF) simulator engines

## Overview

DAF-Sim is a multi-pillar computational framework for modeling and simulating Dissolved Air Flotation (DAF) systems. The simulator integrates fundamental physical models with advanced computational methods to predict and optimize DAF performance in water treatment applications.

## Project Structure

### Components

#### Pillar 2: PBM Equation Shell (Kinetics Solver)

**Status:** ✅ Implemented

**Location:** `pbm_equation_shell/`

**Description:** Discretized Population Balance Model (PBM) for flocculation kinetics based on Abreu et al. (2021). Implements the governing equation for particle size distribution evolution through aggregation and breakage mechanisms.

**Key Features:**
- Complete implementation of discretized PBE (Equation 2)
- Support for aggregation (perikinetic, differential sedimentation, orthokinetic)
- Support for particle breakage
- Comprehensive validation and testing suite
- Mass conservation verification tools

**Documentation:** See `pbm_equation_shell/README.md`

**Quick Start:**
```bash
cd pbm_equation_shell
pip install -r requirements.txt
python tests/test_population_balance_model.py
```

## Installation

```bash
# Clone the repository
git clone <repository-url>
cd DAF-Sim

# Install dependencies for PBM Equation Shell
cd pbm_equation_shell
pip install -r requirements.txt
```

## Usage Example

```python
from pbm_equation_shell import PopulationBalanceModel, construct_beta_matrix
import numpy as np

# Initialize PBM with 5 size classes
pbm = PopulationBalanceModel(n_classes=5)

# Define initial particle distribution
N = np.array([1e15, 5e14, 1e14, 5e13, 1e13])  # #/m³

# Calculate collision frequencies and solve
# See pbm_equation_shell/README.md for complete examples
```

## Development Roadmap

- [x] **Pillar 2:** PBM Equation Shell (Kinetics Solver)
- [ ] **Pillar 1:** Hydrodynamics Module
- [ ] **Pillar 3:** Chemistry and Physics Models
- [ ] **Pillar 4:** Machine Learning Integration
- [ ] **Integration:** Multi-pillar coupling framework

## Testing

Each component includes comprehensive unit tests and integration tests.

```bash
# Test PBM Equation Shell
cd pbm_equation_shell
python -m pytest tests/ -v
```

## Contributing

This is a multi-pillar development project. Each pillar is designed to be modular and independently verifiable while supporting integration with other components.

## License

See LICENSE file for details.

## References

- Abreu, D. A. d. S., et al. (2021). "Integration of first principle models and machine learning in a modeling framework: An application to flocculation" Chemical Engineering Science.
