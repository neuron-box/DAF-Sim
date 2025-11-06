# Pillar 3: Core Physics Model for DAF Simulator

## Overview

This module implements the **collision efficiency factor (α_bp)** calculation for bubble-particle interactions in Dissolved Air Flotation (DAF) systems, based on the mechanistic model by **Han et al. (2001)**.

The model uses **trajectory analysis** combined with **DLVO theory** to predict whether particles will successfully attach to bubbles, accounting for:
- Hydrodynamic forces (drag in bubble flow field)
- van der Waals attractive forces
- Electrostatic double layer repulsive forces
- Gravitational forces

## Scientific Background

### Collision Efficiency Factor (α_bp)

The collision efficiency factor represents the probability that a particle approaching a bubble will successfully attach to it. It is calculated as:

**α_bp = η_collision × α_attachment**

Where:
- **η_collision**: Collision efficiency from trajectory analysis (fraction of particles that physically reach the bubble surface)
- **α_attachment**: Attachment efficiency from DLVO analysis (fraction of colliding particles that overcome the energy barrier to attach)

### Trajectory Analysis Method

The trajectory analysis method solves the equations of motion for a particle in the flow field of a rising bubble:

**m_eff × dv/dt = F_hydrodynamic + F_DLVO + F_gravity**

The particle trajectory is integrated numerically to determine the **critical radial offset** (y_c), which is the maximum initial lateral distance from which a particle can still collide with the bubble.

The collision efficiency is then:

**η_c = (y_c / R_b)²**

### DLVO Theory

The DLVO (Derjaguin-Landau-Verwey-Overbeek) theory describes the total interaction energy between the bubble and particle:

**V_total = V_vdW + V_EDL**

#### 1. Van der Waals Force (Attractive)

For two spheres (particle radius R_p, bubble radius R_b) separated by distance h:

```
F_vdW = -A × R_p × R_b / (6 × h² × (R_p + R_b))

V_vdW = -A × R_p × R_b / (6 × h × (R_p + R_b))
```

Where A is the Hamaker constant [J].

#### 2. Electrostatic Double Layer Force (Repulsive for Like Charges)

```
F_EDL = 2π × ε₀ × εᵣ × κ × R_p × R_b / (R_p + R_b) ×
        [(ζ_p² + ζ_b²) × exp(-κh) - 2ζ_p×ζ_b × exp(-2κh)]

V_EDL = 2π × ε₀ × εᵣ × R_p × R_b / (R_p + R_b) ×
        [(ζ_p² + ζ_b²) × (1/κ) × exp(-κh) + 2ζ_p×ζ_b × (1/2κ) × exp(-2κh)]
```

Where:
- ε₀ = 8.854×10⁻¹² F/m (vacuum permittivity)
- εᵣ = relative dielectric constant (≈78.5 for water)
- κ = Debye parameter (inverse Debye length) [1/m]
- ζ_p, ζ_b = zeta potentials of particle and bubble [V]
- h = surface-to-surface separation [m]

#### 3. Debye Length

The electrical double layer thickness (Debye length) for aqueous solutions at 25°C:

```
λ_D = 0.304 / √I  [nm]
```

Where I is the ionic strength [mol/L].

### Hydrodynamic Forces

The hydrodynamic drag force on the particle follows Stokes law:

```
F_drag = 6π × μ × R_p × v_rel
```

Where:
- μ = dynamic viscosity [Pa·s]
- R_p = particle radius [m]
- v_rel = relative velocity between particle and fluid [m/s]

The flow field around a rising bubble with a mobile interface is approximated by:

```
v_r = U × (1 - (R_b/r)³) × cos(θ)
v_θ = -U × (1 + 0.5 × (R_b/r)³) × sin(θ)
```

Where U is the bubble rise velocity and θ is the angle from the vertical axis.

## Installation

### Requirements

- Python 3.8+
- NumPy
- SciPy

### Install Dependencies

```bash
pip install -r requirements.txt
```

### Install Package

```bash
cd pillar3_physics_model
pip install -e .
```

## Usage

### Basic Example

```python
from pillar3_physics_model import (
    Particle, Bubble, FluidProperties,
    calculate_attachment_efficiency
)

# Define a particle (10 μm diameter, negative zeta potential)
particle = Particle(
    diameter=10e-6,           # [m]
    zeta_potential=-0.025,    # [V]
    density=2650.0,           # [kg/m³]
    hamaker_constant=5e-21    # [J]
)

# Define a bubble (50 μm diameter)
bubble = Bubble(
    diameter=50e-6,           # [m]
    zeta_potential=-0.030     # [V]
)

# Define fluid properties (water at 20°C)
fluid = FluidProperties.water_at_20C(
    ionic_strength=0.001,     # [M]
    ph=7.0
)

# Calculate collision efficiency
result = calculate_attachment_efficiency(particle, bubble, fluid)

# Display results
print(f"Collision efficiency factor: α_bp = {result['alpha_bp']:.4f}")
print(f"Collision efficiency: η_c = {result['eta_collision']:.4f}")
print(f"Attachment efficiency: α = {result['alpha_attachment']:.4f}")
print(f"Energy barrier: {result['energy_barrier_kT']:.1f} kT")
print(f"Critical offset: {result['critical_offset']*1e6:.2f} μm")
print(f"Debye length: {result['debye_length']*1e9:.2f} nm")
```

### Typical DAF Conditions

```python
# After coagulation treatment
particle = Particle(
    diameter=15e-6,            # Flocculated particle
    zeta_potential=-0.010,     # Reduced by coagulant
    density=1200.0,            # Lower density (floc structure)
    hamaker_constant=8e-21     # Hydrophobic after treatment
)

bubble = Bubble(
    diameter=40e-6,            # Typical DAF bubble size
    zeta_potential=-0.020
)

fluid = FluidProperties.water_at_20C(
    ionic_strength=0.01,       # After coagulant addition
    ph=7.0
)

result = calculate_attachment_efficiency(particle, bubble, fluid)
```

## Module Structure

```
pillar3_physics_model/
├── src/
│   └── pillar3_physics_model/
│       ├── __init__.py                 # Package initialization
│       ├── particle.py                 # Particle class
│       ├── bubble.py                   # Bubble class
│       ├── fluid_properties.py         # Fluid properties class
│       ├── dlvo_forces.py              # DLVO force calculations
│       ├── hydrodynamic_forces.py      # Hydrodynamic force calculations
│       ├── trajectory_solver.py        # Trajectory integration
│       └── collision_efficiency.py     # Main calculation function
├── tests/
│   ├── __init__.py
│   └── test_physics_model.py           # Comprehensive unit tests
├── README.md                           # This file
├── requirements.txt                    # Python dependencies
└── setup.py                            # Package setup
```

## Running Tests

Run the comprehensive test suite:

```bash
cd pillar3_physics_model
python -m pytest tests/ -v
```

Or run tests directly:

```bash
python tests/test_physics_model.py
```

## Test Coverage

The test suite includes:

1. **Data Classes**: Validation of Particle, Bubble, and FluidProperties
2. **DLVO Forces**: van der Waals and electrostatic force calculations
3. **Hydrodynamic Forces**: Drag forces and flow field calculations
4. **Trajectory Solver**: Numerical integration and collision detection
5. **Collision Efficiency**: End-to-end calculation validation
6. **Integration Tests**: Realistic DAF scenarios

## Physical Constants

The following physical constants are used:

- Vacuum permittivity: ε₀ = 8.854×10⁻¹² F/m
- Boltzmann constant: k_B = 1.381×10⁻²³ J/K
- Elementary charge: e = 1.602×10⁻¹⁹ C
- Gravitational acceleration: g = 9.81 m/s²

## Typical Parameter Ranges

### Particles
- Diameter: 1-100 μm
- Zeta potential: -50 to +50 mV
- Density: 1000-3000 kg/m³
- Hamaker constant: 0.5-10 × 10⁻²¹ J

### Bubbles
- Diameter: 10-100 μm (DAF typical)
- Zeta potential: -20 to -50 mV (for air-water interface)

### Fluid
- Temperature: 278-298 K (5-25°C)
- Ionic strength: 0.0001-0.1 M
- pH: 4-10

## References

1. **Han, M.Y., Kim, W., & Dockko, S.** (2001). "Collision efficiency factor of bubble and particle (α_bp) in DAF: theory and experimental verification." *Water Science and Technology*, 43(8), 139-144.

2. **Han, M., & Lawler, D. F.** (1991). "The (Relative) Insignificance of G in Flocculation." *Journal AWWA*, 83(3), 71-77.

3. **Israelachvili, J. N.** (2011). *Intermolecular and Surface Forces* (3rd ed.). Academic Press.

4. **Yoon, R. H., & Luttrell, G. H.** (1989). "The Effect of Bubble Size on Fine Particle Flotation." *Mineral Processing and Extractive Metallurgy Review*, 5(1-4), 101-122.

5. **Edzwald, J. K.** (2010). "Dissolved air flotation and me." *Water Research*, 44(7), 2077-2106.

## Author

Developed for the DAF-Sim multi-pillar simulator project.

## License

This code is provided as-is for scientific and educational purposes.
