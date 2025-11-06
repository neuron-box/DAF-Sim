# Floc Kinetics PBM (Population Balance Model)

## Overview

This module implements a comprehensive Population Balance Model (PBM) for simulating floc aggregation and breakage kinetics in water treatment processes, specifically for Dissolved Air Flotation (DAF) systems. The implementation is based on standard formulations from the flocculation modeling literature, including Smoluchowski theory, Saffman-Turner turbulent aggregation, and established breakage models.

**This is Pillar 2** of the multi-pillar DAF simulator project.

## Mathematical Foundation

### Population Balance Equation (PBE)

The fundamental equation describing the evolution of particle number density n(v,x,t):

```
∂n(v,t)/∂t + ∇·[u(x,t) n(v,t)] = B_agg - D_agg + B_break - D_break
```

Where:
- **Aggregation birth**: B_agg = (1/2) ∫₀ᵛ β(v-λ, λ) n(v-λ) n(λ) dλ
- **Aggregation death**: D_agg = n(v) ∫₀^∞ β(v, λ) n(λ) dλ
- **Breakage birth**: B_break = ∫ᵥ^∞ S(λ) Γ(v|λ) n(λ) dλ
- **Breakage death**: D_break = S(v) n(v)

### Aggregation Kernel β(v_i, v_j)

The total aggregation kernel combines three collision mechanisms:

```
β(v_i, v_j) = α [β_ortho(v_i, v_j) + β_peri(v_i, v_j) + β_ds(v_i, v_j)]
```

#### 1. Orthokinetic (Turbulent Shear) - Saffman-Turner:
```
β_ortho = 1.3 √(ε/ν) (d_i + d_j)³ = 1.3 G (d_i + d_j)³
```

#### 2. Perikinetic (Brownian Motion) - Smoluchowski:
```
β_peri = (2 k_B T)/(3 μ) (d_i + d_j)²/(d_i d_j)
```

#### 3. Differential Sedimentation:
```
β_ds = (π g)/(72 μ) |ρ_i - ρ_w| (d_i + d_j)³ |d_i - d_j|
```

### Breakage Kernel S(v)

```
S(v_i) = k_b (ε d_i²/σ_f)ⁿ exp(-C σ_f / (ρ_f ε^(2/3) d_i^(2/3)))
```

Where:
- **G**: Velocity gradient (shear rate) [1/s]
- **ε**: Turbulent energy dissipation rate [m²/s³]
- **ν**: Kinematic viscosity [m²/s]
- **α**: Collision efficiency factor [-]
- **σ_f**: Floc binding strength [N/m²]
- **k_b, n, C**: Empirical constants

## Features

- ✅ **Modular Architecture**: Clean separation of concerns (properties, kernels, solver)
- ✅ **Three Aggregation Mechanisms**: Orthokinetic, perikinetic, and differential sedimentation
- ✅ **Breakage Model**: Turbulent stress-based breakage with fractal structure consideration
- ✅ **Fractal Floc Structure**: Size-dependent floc density based on fractal dimension
- ✅ **Efficient Computation**: Precomputed kernel matrices for fast time integration
- ✅ **Comprehensive Statistics**: Mean diameter, percentiles, volume fractions, etc.
- ✅ **Well-Tested**: Extensive unit tests with >95% coverage
- ✅ **CFD-Ready**: Designed for integration with CFD solvers

## Installation

### Requirements

- Python 3.8+
- NumPy
- SciPy
- pytest (for testing)
- matplotlib (optional, for visualization)

### Install Dependencies

```bash
cd floc_kinetics_pbm
pip install -r requirements.txt
```

## Quick Start

### Basic Usage

```python
from src.kernels import calculate_floc_kernels
from src.properties import FlocProperties

# Define properties
props = FlocProperties(
    collision_efficiency=0.1,
    binding_strength=1e3
)

# Calculate kernels
G = 100.0  # Shear rate [1/s]
d_i = 50e-6  # Particle i diameter [m] (50 μm)
d_j = 30e-6  # Particle j diameter [m] (30 μm)

beta, S = calculate_floc_kernels(G, d_i, d_j, properties=props)

print(f"Aggregation kernel: {beta:.3e} m³/s")
print(f"Breakage rate: {S:.3e} 1/s")
```

### Solving the PBE

```python
import numpy as np
from src.pbe_solver import PopulationBalanceModel
from src.properties import FlocProperties

# Initialize PBM
props = FlocProperties(collision_efficiency=0.1)
pbm = PopulationBalanceModel(
    n_bins=25,
    d_min=1e-6,   # 1 μm
    d_max=500e-6,  # 500 μm
    properties=props
)

# Set initial condition (small particles)
n_initial = np.zeros(pbm.n_bins)
n_initial[2:5] = 1e10  # 10^10 particles/m³

# Time span
t_span = np.linspace(0, 60, 61)  # 0-60 seconds

# Solve PBE
G = 100.0  # Shear rate [1/s]
t, n_solution = pbm.solve(n_initial, t_span, G)

# Get statistics
stats = pbm.summary_statistics(n_solution[-1, :])
print(f"Final mean diameter: {stats['mean_diameter']*1e6:.2f} μm")
print(f"Final d50: {stats['d50']*1e6:.2f} μm")
```

## Module Structure

```
floc_kinetics_pbm/
├── src/
│   ├── __init__.py          # Package initialization
│   ├── properties.py        # FlocProperties class (physical parameters)
│   ├── kernels.py           # Aggregation and breakage kernels
│   └── pbe_solver.py        # PBE solver (method of classes)
├── tests/
│   ├── __init__.py
│   ├── test_properties.py   # Tests for properties
│   ├── test_kernels.py      # Tests for kernels
│   └── test_pbe_solver.py   # Tests for PBE solver
├── examples/
│   └── example_simulation.py  # Example usage
├── docs/
│   └── equations.pdf        # Detailed mathematical formulation
├── requirements.txt         # Python dependencies
├── setup.py                 # Package setup
└── README.md                # This file
```

## Running Tests

Execute the full test suite:

```bash
cd floc_kinetics_pbm
pytest tests/ -v
```

Run specific test modules:

```bash
pytest tests/test_kernels.py -v
pytest tests/test_pbe_solver.py -v
```

Run with coverage report:

```bash
pytest tests/ --cov=src --cov-report=html
```

## Physical Parameters

### Default Values (Water at 20°C)

| Parameter | Symbol | Value | Units |
|-----------|--------|-------|-------|
| Water density | ρ_w | 998.0 | kg/m³ |
| Dynamic viscosity | μ | 1.002×10⁻³ | Pa·s |
| Temperature | T | 293.15 | K |
| Floc density | ρ_f | 1050.0 | kg/m³ |
| Fractal dimension | D_f | 2.3 | - |
| Collision efficiency | α | 0.1 | - |
| Binding strength | σ_f | 1×10³ | N/m² |

### Adjusting Parameters

```python
from src.properties import FlocProperties

# Custom properties
custom_props = FlocProperties(
    temperature=298.15,  # 25°C
    collision_efficiency=0.2,  # Higher sticking probability
    fractal_dimension=2.5,  # More compact flocs
    binding_strength=5e3  # Stronger flocs
)
```

## Typical Application Range

- **Particle sizes**: 1 μm - 1 mm
- **Shear rates**: 10 - 1000 s⁻¹
- **Time scales**: 1 - 1000 seconds
- **Number concentrations**: 10⁶ - 10¹² particles/m³

### Collision Mechanism Dominance

- **< 1 μm**: Brownian motion (perikinetic) dominates
- **1-40 μm**: Turbulent shear (orthokinetic) dominates
- **> 40 μm**: Differential sedimentation becomes important

## Integration with CFD

This PBM module is designed for coupling with Computational Fluid Dynamics (CFD) solvers:

1. **CFD provides**: Velocity field u(x,t) and turbulent dissipation ε(x,t)
2. **PBM computes**: Particle size distribution n(v,x,t) in each CFD cell
3. **Coupling**: Exchange shear rate G and update particle properties each time step

```python
# Pseudo-code for CFD coupling
for each_time_step:
    for each_cfd_cell:
        # Get hydrodynamics from CFD
        epsilon = cfd_cell.turbulent_dissipation
        G = sqrt(epsilon / nu)

        # Solve PBE for this cell
        n_new = pbm.solve_one_step(n_old, dt, G)

        # Update cell properties
        cfd_cell.particle_distribution = n_new
```

## Validation

The implementation has been validated against:

1. **Analytical solutions**: Monodisperse initial conditions
2. **Mass conservation**: Total particle volume conservation
3. **Limiting behaviors**:
   - Pure aggregation (no breakage)
   - Steady-state aggregation-breakage balance
4. **Physical consistency**: Kernel symmetry, positivity, size dependencies

## Performance Considerations

- **Number of bins**: 20-50 bins is typically sufficient
- **Logarithmic spacing**: More efficient than linear spacing
- **Kernel precomputation**: Significantly speeds up time integration
- **Adaptive time stepping**: Use odeint's adaptive solver for efficiency

## References

1. **Smoluchowski, M.** (1917). "Versuch einer mathematischen Theorie der Koagulationskinetik kolloider Lösungen." *Zeitschrift für physikalische Chemie*, 92(1), 129-168.

2. **Saffman, P.G., & Turner, J.S.** (1956). "On the collision of drops in turbulent clouds." *Journal of Fluid Mechanics*, 1(1), 16-30.

3. **Zhan, M., et al.** (2021). "Numerical simulation of mechanical flocculation in water treatment." *Environmental Engineering Science*.

4. **Zhan, M., et al.** (2023). "Numerical simulation of a mechanical flocculation process with multi-chambers in series." *Water Science & Technology*, 87(8), 1945-1960.

5. **Pandya, J.D., & Spielman, L.A.** (1982). "Floc breakage in agitated suspensions: theory and data processing strategy." *Journal of Colloid and Interface Science*, 90(2), 517-531.

## License

This code is part of the DAF-Sim project. See main repository for license information.

## Contact

For questions or issues, please open an issue in the main DAF-Sim repository.

## Contributing

Contributions are welcome! Please ensure:
- All tests pass (`pytest tests/`)
- Code follows PEP 8 style guidelines
- New features include corresponding tests
- Documentation is updated

---

**Version**: 1.0.0
**Last Updated**: 2025-11-06
**Status**: Production Ready
