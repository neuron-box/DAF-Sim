# PBM Equation Shell - Population Balance Model for Flocculation

## Overview

This component implements the discretized **Population Balance Equation (PBE)** for modeling particle flocculation in Dissolved Air Flotation (DAF) systems. The implementation is based on the explicit, discretized PBM framework described in:

> **Abreu, D. A. d. S., et al. (2021).** "Integration of first principle models and machine learning in a modeling framework: An application to flocculation" *Chemical Engineering Science*.

This module serves as the core "equation shell" for the DAF simulator's Pillar 2 kinetics solver, providing the fundamental mathematical framework for tracking particle size distribution evolution through aggregation and breakage processes.

## Mathematical Foundation

### The Discretized Population Balance Equation

The discretized PBE (Equation 2 from Abreu et al., 2021) describes the time evolution of particle number concentrations in each size class:

```
dNᵢ/dt = (Birth by aggregation) - (Death by aggregation)
         + (Birth by breakage) - (Death by breakage)
```

**Explicitly:**

```
dNᵢ/dt = (1/2) · Σⱼ₌₁ⁱ⁻¹ βⱼ,ᵢ₋ⱼ · αⱼ,ᵢ₋ⱼ · Nⱼ · Nᵢ₋ⱼ
         - Nᵢ · Σⱼ₌₁ⁿ βᵢ,ⱼ · αᵢ,ⱼ · Nⱼ
         + Σⱼ₌ᵢ₊₁ⁿ Sⱼ · γⱼ,ᵢ · Nⱼ
         - Sᵢ · Nᵢ
```

### Term Definitions

#### State Variables
- **Nᵢ** : Number concentration of particles in size class *i* [#/m³]

#### Aggregation Kernels
- **βᵢ,ⱼ** : Collision frequency between particles of size class *i* and *j* [m³/s]
  - Constructed from three physical mechanisms:
    1. **Perikinetic** (Brownian motion): Important for small particles (< 1 μm)
    2. **Differential sedimentation**: Dominant for larger particles with density differences
    3. **Orthokinetic** (shear-induced): Significant in turbulent flows

- **αᵢ,ⱼ** : Collision efficiency between particles of size class *i* and *j* [-]
  - Dimensionless probability that a collision results in permanent attachment
  - Ranges from 0 (no attachment) to 1 (perfect attachment)
  - Accounts for surface chemistry, repulsive forces, etc.

#### Breakage Kernels
- **Sᵢ** : Breakage rate of particles in size class *i* [1/s]
  - Frequency at which particles of size *i* break apart

- **γⱼ,ᵢ** : Breakage distribution function [-]
  - Probability that breakage of particle *j* produces particle *i*
  - Satisfies normalization: Σᵢ γⱼ,ᵢ = 1

### Collision Frequency Construction

The collision frequency β is computed as the sum of three mechanisms:

**Perikinetic (Brownian motion):**
```
βₚ = (2kᵦT)/(3μ) · (dᵢ + dⱼ)²/(dᵢ · dⱼ)
```

**Differential sedimentation:**
```
βDS = (πg|Δρ|)/(72μ) · (dᵢ + dⱼ)³ · |dᵢ - dⱼ|
```

**Orthokinetic (shear):**
```
βO = (4/3) · G · (dᵢ + dⱼ)³
```

Where:
- kᵦ = Boltzmann constant [1.381×10⁻²³ J/K]
- T = Temperature [K]
- μ = Dynamic viscosity [Pa·s]
- g = Gravitational acceleration [9.81 m/s²]
- Δρ = Density difference [kg/m³]
- G = Shear rate [1/s]
- dᵢ, dⱼ = Particle diameters [m]

## Installation

### Requirements

- Python ≥ 3.8
- NumPy ≥ 1.20

### Setup

```bash
# Clone the repository
cd DAF-Sim

# Install dependencies
pip install numpy

# For development (testing)
pip install pytest
```

## Usage

### Basic Example

```python
import numpy as np
from pbm_equation_shell import PopulationBalanceModel, construct_beta_matrix

# Define system parameters
n_classes = 5  # Number of size classes

# Initialize the PBM
pbm = PopulationBalanceModel(n_classes=n_classes)

# Initial particle distribution (mostly small particles)
N = np.array([1e15, 5e14, 1e14, 5e13, 1e13])  # #/m³

# Define particle sizes (1 to 100 μm)
particle_diameters = np.logspace(-6, -4, n_classes)  # m

# Construct collision frequency matrix
beta_matrix = construct_beta_matrix(
    n_classes,
    particle_diameters,
    temperature=298.15,      # K
    viscosity=8.9e-4,        # Pa·s (water at 25°C)
    shear_rate=50.0          # 1/s
)

# Define collision efficiency (30% attachment probability)
alpha_matrix = np.ones((n_classes, n_classes)) * 0.3

# No breakage for this example
S_vector = np.zeros(n_classes)
gamma_matrix = np.zeros((n_classes, n_classes))

# Calculate time derivatives
dndt = pbm.calculate_dndt(N, alpha_matrix, beta_matrix, S_vector, gamma_matrix)

print("Time derivatives dN/dt:")
for i, rate in enumerate(dndt):
    print(f"  Class {i}: {rate:.2e} #/(m³·s)")
```

### With Breakage

```python
# Enable breakage for larger particles
S_vector = np.array([0.0, 0.0, 0.1, 0.5, 1.0])  # 1/s

# Breakage produces smaller particles
gamma_matrix = np.zeros((n_classes, n_classes))
gamma_matrix[2, 1] = 1.0  # Class 2 breaks to class 1
gamma_matrix[3, 1] = 0.5  # Class 3 breaks to class 1 (50%)
gamma_matrix[3, 2] = 0.5  # Class 3 breaks to class 2 (50%)
gamma_matrix[4, 2] = 0.7  # Class 4 breaks to class 2 (70%)
gamma_matrix[4, 3] = 0.3  # Class 4 breaks to class 3 (30%)

# Calculate with both aggregation and breakage
dndt = pbm.calculate_dndt(N, alpha_matrix, beta_matrix, S_vector, gamma_matrix)
```

### Mass Conservation Check

```python
# Calculate particle volumes
particle_volumes = (4/3) * np.pi * (particle_diameters/2)**3

# Validate mass conservation
is_conserved, mass_rate = pbm.validate_mass_conservation(
    N, dndt, particle_volumes, tolerance=1e-6
)

print(f"Mass conserved: {is_conserved}")
print(f"Mass rate: {mass_rate:.2e} m³/(m³·s)")
```

## Running Tests

The module includes comprehensive unit tests and integration tests:

```bash
# Run all tests
cd pbm_equation_shell
python -m pytest tests/ -v

# Or run directly
python tests/test_population_balance_model.py
```

### Test Coverage

The test suite includes:

1. **Initialization tests**: Verify correct setup and validation
2. **Shape tests**: Ensure correct array dimensions
3. **Physics tests**:
   - Zero aggregation/breakage → zero derivatives
   - Pure aggregation behavior
   - Pure breakage behavior
   - Aggregation-breakage balance
4. **Input validation tests**: All kernels and vectors
5. **Symmetry tests**: Kernel symmetry properties
6. **Mass conservation tests**: Physical consistency
7. **Integration tests**: Realistic flocculation scenarios

Expected output:
```
======================================================================
SIMPLE FLOCCULATION EXAMPLE
======================================================================

Initial particle concentrations (#/m³):
  Class 0: 1.00e+15
  Class 1: 5.00e+14
  ...

RUNNING UNIT TESTS
======================================================================
test_aggregation_breakage_balance ... ok
test_beta_matrix_positive ... ok
...
Ran 20 tests in 0.XXXs

OK
```

## API Reference

### `PopulationBalanceModel`

**Constructor:**
```python
PopulationBalanceModel(n_classes: int, validate_inputs: bool = True)
```

**Methods:**

#### `calculate_dndt(N, alpha_matrix, beta_matrix, S_vector, gamma_matrix)`
Calculate time derivatives dN/dt for all size classes.

**Parameters:**
- `N` (ndarray): Particle concentrations [#/m³], shape `(n_classes,)`
- `alpha_matrix` (ndarray): Collision efficiencies [-], shape `(n_classes, n_classes)`
- `beta_matrix` (ndarray): Collision frequencies [m³/s], shape `(n_classes, n_classes)`
- `S_vector` (ndarray): Breakage rates [1/s], shape `(n_classes,)`
- `gamma_matrix` (ndarray): Breakage distributions [-], shape `(n_classes, n_classes)`

**Returns:**
- `dndt` (ndarray): Time derivatives [#/(m³·s)], shape `(n_classes,)`

#### `validate_mass_conservation(N, dndt, particle_volumes, tolerance=1e-10)`
Check mass conservation in the system.

**Parameters:**
- `N` (ndarray): Current concentrations [#/m³]
- `dndt` (ndarray): Time derivatives [#/(m³·s)]
- `particle_volumes` (ndarray): Particle volumes [m³], shape `(n_classes,)`
- `tolerance` (float): Conservation tolerance

**Returns:**
- `is_conserved` (bool): Whether mass is conserved
- `mass_rate` (float): Rate of mass change [m³/(m³·s)]

### `construct_beta_matrix()`

**Function:**
```python
construct_beta_matrix(
    n_classes,
    particle_diameters,
    temperature=298.15,
    viscosity=8.9e-4,
    density_fluid=998.2,
    density_particle=1050.0,
    shear_rate=0.0
)
```

Construct collision frequency matrix from physical mechanisms.

**Parameters:**
- `n_classes` (int): Number of size classes
- `particle_diameters` (ndarray): Diameters [m], shape `(n_classes,)`
- `temperature` (float): Temperature [K], default 298.15
- `viscosity` (float): Viscosity [Pa·s], default 8.9e-4
- `density_fluid` (float): Fluid density [kg/m³], default 998.2
- `density_particle` (float): Particle density [kg/m³], default 1050.0
- `shear_rate` (float): Shear rate [1/s], default 0.0

**Returns:**
- `beta_matrix` (ndarray): Collision frequencies [m³/s], shape `(n_classes, n_classes)`

## Integration with DAF Simulator

This PBM equation shell is designed to integrate with the larger DAF simulator as **Pillar 2: Kinetics Solver**. It provides:

1. **Modular interface**: Clean separation between kinetics (this module) and other pillars
2. **Flexible kernels**: Users provide α, β, S, γ to customize physics
3. **Validation tools**: Built-in mass conservation and input checking
4. **Performance**: Optimized NumPy operations for speed

### Integration Pattern

```python
# In DAF simulator main loop
from pbm_equation_shell import PopulationBalanceModel

# Initialize once
pbm = PopulationBalanceModel(n_classes=100)

# In time-stepping loop
for t in time_steps:
    # Get current state from Pillar 1 (hydrodynamics)
    N_current = get_particle_distribution()

    # Get kernels from Pillar 3 (chemistry/physics models)
    alpha = calculate_collision_efficiency(...)
    beta = calculate_collision_frequency(...)
    S, gamma = calculate_breakage_kernels(...)

    # Compute derivatives (this module)
    dndt = pbm.calculate_dndt(N_current, alpha, beta, S, gamma)

    # Update state
    N_next = time_integrator(N_current, dndt, dt)
```

## Design Decisions

### Why This Structure?

1. **Explicit summations**: The implementation uses explicit loops for clarity and correspondence to the mathematical equation. This makes the code directly verifiable against Equation 2 in the paper.

2. **Matrix inputs**: Kernels are provided as matrices/vectors rather than functions, allowing:
   - Pre-computation of expensive operations
   - Easier integration with ML models (Pillar 4)
   - Clear separation between physics models and kinetics solver

3. **Validation toggles**: Input validation can be disabled for performance in production runs after debugging.

4. **Zero-indexed arrays**: Python uses 0-indexing, so size class *i* in the equations corresponds to `N[i]` in code.

### Performance Considerations

For large systems (n_classes > 100), consider:
- Sparse matrix representations for kernels
- Vectorization of aggregation terms
- JIT compilation (e.g., Numba)
- GPU acceleration for matrix operations

The current implementation prioritizes clarity and correctness over raw speed.

## Limitations and Assumptions

1. **Fixed size classes**: Number and bounds of size classes must be defined a priori
2. **Symmetric kernels**: `construct_beta_matrix()` assumes symmetric collision kernels
3. **Well-mixed**: Assumes spatially homogeneous distribution (appropriate for CSTR)
4. **No external forces**: Beyond gravity in differential sedimentation
5. **Binary aggregation**: No triple collisions or higher-order events

## Citation

If you use this code in research, please cite:

```bibtex
@article{abreu2021integration,
  title={Integration of first principle models and machine learning in a modeling framework: An application to flocculation},
  author={Abreu, D. A. d. S. and others},
  journal={Chemical Engineering Science},
  year={2021}
}
```

## License

Part of the DAF-Sim project. See main repository for license details.

## Contributing

This module is part of the DAF-Sim multi-pillar architecture. For changes:
1. Ensure all tests pass
2. Maintain mathematical correspondence to Equation 2 (Abreu et al., 2021)
3. Add tests for new features
4. Update this README for API changes

## Contact

For questions or issues specific to this component, please open an issue in the main DAF-Sim repository.

---

**Version:** 1.0.0
**Last Updated:** November 2025
**Component:** Pillar 2 - Kinetics Solver (PBM Equation Shell)
