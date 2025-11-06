# Discretized Population Balance Equation

## Complete Mathematical Formulation

This document provides the complete discretized Population Balance Equation (PBE) as extracted from Abreu et al. (2021), Equation 2.

### The Discretized PBE (Equation 2)

For each size class $i$ where $i = 1, 2, \ldots, n$, the time evolution of particle number concentration is governed by:

$$
\frac{dN_i}{dt} = \frac{1}{2} \sum_{j=1}^{i-1} \beta_{j,i-j} \alpha_{j,i-j} N_j N_{i-j} - N_i \sum_{j=1}^{n} \beta_{i,j} \alpha_{i,j} N_j + \sum_{j=i+1}^{n} S_j \gamma_{j,i} N_j - S_i N_i
$$

### Term-by-Term Breakdown

#### Term 1: Birth by Aggregation

$$
\frac{1}{2} \sum_{j=1}^{i-1} \beta_{j,i-j} \alpha_{j,i-j} N_j N_{i-j}
$$

**Physical meaning:** Particles of size class $i$ are formed by the collision and aggregation of two smaller particles of sizes $j$ and $i-j$.

**Components:**
- $\beta_{j,i-j}$: Collision frequency between particles of size $j$ and $i-j$ [m³/s]
- $\alpha_{j,i-j}$: Collision efficiency (attachment probability) [-]
- $N_j$, $N_{i-j}$: Number concentrations of reactant particles [#/m³]
- Factor $\frac{1}{2}$: Accounts for double-counting (since $j + (i-j) = (i-j) + j$)

#### Term 2: Death by Aggregation

$$
-N_i \sum_{j=1}^{n} \beta_{i,j} \alpha_{i,j} N_j
$$

**Physical meaning:** Particles of size class $i$ are lost through collision with any other particle $j$, forming a larger aggregate.

**Components:**
- $\beta_{i,j}$: Collision frequency between particles of size $i$ and $j$ [m³/s]
- $\alpha_{i,j}$: Collision efficiency [-]
- $N_i$: Number concentration of size $i$ (being depleted) [#/m³]
- $N_j$: Number concentration of collision partner [#/m³]

#### Term 3: Birth by Breakage

$$
\sum_{j=i+1}^{n} S_j \gamma_{j,i} N_j
$$

**Physical meaning:** Particles of size class $i$ are formed by the breakage of larger particles $j > i$.

**Components:**
- $S_j$: Breakage rate of particles in size class $j$ [1/s]
- $\gamma_{j,i}$: Breakage distribution function (probability that breakage of $j$ produces $i$) [-]
- $N_j$: Number concentration of parent particles [#/m³]
- Constraint: $\sum_{i=1}^{j} \gamma_{j,i} = 1$ (all fragments are accounted for)

#### Term 4: Death by Breakage

$$
-S_i N_i
$$

**Physical meaning:** Particles of size class $i$ are lost through their own breakage into smaller fragments.

**Components:**
- $S_i$: Breakage rate of size class $i$ [1/s]
- $N_i$: Number concentration being depleted [#/m³]

## Variable Definitions

### State Variable

| Symbol | Description | Units |
|--------|-------------|-------|
| $N_i$ | Number concentration of particles in size class $i$ | #/m³ |

### Aggregation Kernels

| Symbol | Description | Units | Range |
|--------|-------------|-------|-------|
| $\beta_{i,j}$ | Collision frequency between size classes $i$ and $j$ | m³/s | $[0, \infty)$ |
| $\alpha_{i,j}$ | Collision efficiency between size classes $i$ and $j$ | - | $[0, 1]$ |

### Breakage Kernels

| Symbol | Description | Units | Range |
|--------|-------------|-------|-------|
| $S_i$ | Breakage rate of particles in size class $i$ | 1/s | $[0, \infty)$ |
| $\gamma_{j,i}$ | Probability that breakage of $j$ produces $i$ | - | $[0, 1]$ |

## Collision Frequency Function $\beta_{i,j}$

The collision frequency is constructed from three physical mechanisms:

$$
\beta_{i,j} = \beta_{i,j}^{\text{peri}} + \beta_{i,j}^{\text{sed}} + \beta_{i,j}^{\text{ortho}}
$$

### Perikinetic (Brownian Motion)

$$
\beta_{i,j}^{\text{peri}} = \frac{2 k_B T}{3 \mu} \frac{(d_i + d_j)^2}{d_i d_j}
$$

**Dominant for:** Small particles ($d < 1$ μm) in quiescent or low-shear conditions

**Parameters:**
- $k_B = 1.380649 \times 10^{-23}$ J/K (Boltzmann constant)
- $T$: Absolute temperature [K]
- $\mu$: Dynamic viscosity of fluid [Pa·s]
- $d_i$, $d_j$: Particle diameters [m]

### Differential Sedimentation

$$
\beta_{i,j}^{\text{sed}} = \frac{\pi g |\Delta \rho|}{72 \mu} (d_i + d_j)^3 |d_i - d_j|
$$

**Dominant for:** Larger particles with different sizes settling at different velocities

**Parameters:**
- $g = 9.81$ m/s² (gravitational acceleration)
- $\Delta \rho = |\rho_p - \rho_f|$: Density difference between particle and fluid [kg/m³]
- $\mu$: Dynamic viscosity [Pa·s]
- $d_i$, $d_j$: Particle diameters [m]

**Note:** Proportional to $(d_i + d_j)^3 |d_i - d_j|$, so particles of very different sizes collide most frequently.

### Orthokinetic (Shear-Induced)

$$
\beta_{i,j}^{\text{ortho}} = \frac{4}{3} G (d_i + d_j)^3
$$

**Dominant for:** Turbulent or shear flows, typical in flocculation reactors

**Parameters:**
- $G$: Shear rate [1/s]
- $d_i$, $d_j$: Particle diameters [m]

**Note:** Proportional to $(d_i + d_j)^3$, so larger particles collide much more frequently.

## Conservation Properties

### Particle Number Conservation (Pure Aggregation)

For a system with only aggregation ($S_i = 0$ for all $i$), the total number of particles decreases:

$$
\frac{d}{dt} \sum_{i=1}^{n} N_i \leq 0
$$

### Mass Conservation (Pure Aggregation)

For pure aggregation, the total mass is conserved:

$$
\frac{d}{dt} \sum_{i=1}^{n} N_i v_i = 0
$$

where $v_i$ is the volume of particles in size class $i$.

**Proof:** Each aggregation event removes two particles and creates one, conserving total volume.

### Breakage Distribution Constraint

The breakage distribution must satisfy:

$$
\sum_{i=1}^{j} \gamma_{j,i} = 1 \quad \forall j
$$

This ensures all fragments from breakage of particle $j$ are accounted for.

## Dimensionless Analysis

### Typical Scales

For water at 25°C treating flocs in a DAF system:

| Parameter | Typical Value | Units |
|-----------|---------------|-------|
| Temperature $T$ | 298 | K |
| Viscosity $\mu$ | $8.9 \times 10^{-4}$ | Pa·s |
| Density (water) $\rho_f$ | 998 | kg/m³ |
| Density (flocs) $\rho_p$ | 1050 | kg/m³ |
| Shear rate $G$ | 10-100 | 1/s |
| Particle diameter $d$ | $10^{-6}$ to $10^{-3}$ | m |

### Characteristic Times

**Aggregation time scale:**
$$
\tau_{\text{agg}} \sim \frac{1}{\beta \alpha N}
$$

**Breakage time scale:**
$$
\tau_{\text{break}} \sim \frac{1}{S}
$$

**Damköhler number (aggregation vs breakage):**
$$
Da = \frac{\tau_{\text{break}}}{\tau_{\text{agg}}} = \frac{\beta \alpha N}{S}
$$

- $Da \gg 1$: Aggregation-dominated
- $Da \ll 1$: Breakage-dominated
- $Da \sim 1$: Dynamic equilibrium between aggregation and breakage

## Numerical Implementation Notes

### Index Conventions

In the mathematical formulation, indices run from $i = 1$ to $n$.

In the Python implementation, arrays are zero-indexed: $i = 0$ to $n-1$.

**Translation:**
- Math: $N_1, N_2, \ldots, N_n$
- Code: `N[0], N[1], ..., N[n-1]`
- Math: $\sum_{j=1}^{i-1}$
- Code: `for j in range(i):`

### Computational Complexity

For $n$ size classes, the computational cost per time step is:

- **Aggregation terms:** $\mathcal{O}(n^2)$ operations
- **Breakage terms:** $\mathcal{O}(n^2)$ operations
- **Total:** $\mathcal{O}(n^2)$ per time step

For large systems, consider:
- Sparse representations if many kernel values are zero
- Vectorization to leverage BLAS/LAPACK
- Parallel computation for independent size classes

## References

1. **Abreu, D. A. d. S., et al.** (2021). "Integration of first principle models and machine learning in a modeling framework: An application to flocculation." *Chemical Engineering Science*.

2. **Smoluchowski, M. v.** (1917). "Versuch einer mathematischen Theorie der Koagulationskinetik kolloider Lösungen." *Zeitschrift für physikalische Chemie*, 92(1), 129-168.

3. **Saffman, P. G., & Turner, J. S.** (1956). "On the collision of drops in turbulent clouds." *Journal of Fluid Mechanics*, 1(1), 16-30.

4. **Ramkrishna, D.** (2000). *Population Balances: Theory and Applications to Particulate Systems in Engineering.* Academic Press.

---

**Document Version:** 1.0
**Last Updated:** November 2025
**Corresponds to Code Version:** 1.0.0
