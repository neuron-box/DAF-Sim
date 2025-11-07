# Mathematical Formulations for Floc Kinetics PBM

## Population Balance Equation (PBE)

### Complete Form

The Population Balance Equation describes the evolution of the particle number density $n(v,\mathbf{x},t)$ where:
- $v$: particle volume [m³]
- $\mathbf{x}$: spatial position [m]
- $t$: time [s]

$$\frac{\partial n(v,t)}{\partial t} + \nabla \cdot [\mathbf{u}(\mathbf{x},t) n(v,t)] = B_{agg} - D_{agg} + B_{break} - D_{break}$$

### Aggregation Terms

**Birth by Aggregation** (particles of size $v$ formed by collision of smaller particles):

$$B_{agg}(v) = \frac{1}{2} \int_0^v \beta(v-\lambda, \lambda) \, n(v-\lambda) \, n(\lambda) \, d\lambda$$

**Death by Aggregation** (particles of size $v$ collide with other particles):

$$D_{agg}(v) = n(v) \int_0^{\infty} \beta(v, \lambda) \, n(\lambda) \, d\lambda$$

### Breakage Terms

**Birth by Breakage** (particles of size $v$ formed from breakage of larger particles):

$$B_{break}(v) = \int_v^{\infty} S(\lambda) \, \Gamma(v|\lambda) \, n(\lambda) \, d\lambda$$

**Death by Breakage** (particles of size $v$ break into smaller fragments):

$$D_{break}(v) = S(v) \, n(v)$$

### Nomenclature

- $\beta(v_i, v_j)$: Aggregation kernel [m³/s]
- $S(v)$: Breakage rate [1/s]
- $\Gamma(v|v_{parent})$: Daughter particle distribution function [1/m³]

---

## Aggregation Kernel $\beta(v_i, v_j)$

The total aggregation kernel combines three collision mechanisms:

$$\beta(v_i, v_j) = \alpha \left[\beta_{ortho}(v_i, v_j) + \beta_{peri}(v_i, v_j) + \beta_{ds}(v_i, v_j)\right]$$

where $\alpha$ is the collision efficiency (sticking probability), $0 \leq \alpha \leq 1$.

### 1. Orthokinetic Aggregation (Turbulent Shear)

Based on the **Saffman-Turner** formulation for turbulent collision in homogeneous isotropic turbulence:

$$\beta_{ortho}(v_i, v_j) = 1.3 \sqrt{\frac{\varepsilon}{\nu}} (d_i + d_j)^3 = 1.3 \, G \, (d_i + d_j)^3$$

where:
- $\varepsilon$: Turbulent energy dissipation rate [m²/s³]
- $\nu$: Kinematic viscosity [m²/s]
- $G = \sqrt{\varepsilon/\nu}$: Shear rate (velocity gradient) [1/s]
- $d_i, d_j$: Equivalent spherical diameters of particles $i$ and $j$ [m]

**Relation to volume:**
$$d = \left(\frac{6v}{\pi}\right)^{1/3}$$

**Physical interpretation:** In turbulent flow, particles are brought together by velocity gradients (eddies). The collision frequency scales with $G$ and the cubic power of the collision cross-section $(d_i + d_j)^3$.

**Dominant regime:** Particles in the size range 1-40 μm in turbulent flow.

### 2. Perikinetic Aggregation (Brownian Motion)

Based on **Smoluchowski's** theory for Brownian collisions:

$$\beta_{peri}(v_i, v_j) = \frac{2 k_B T}{3 \mu} \frac{(d_i + d_j)^2}{d_i \, d_j}$$

where:
- $k_B = 1.380649 \times 10^{-23}$ J/K: Boltzmann constant
- $T$: Absolute temperature [K]
- $\mu$: Dynamic viscosity [Pa·s]

**Physical interpretation:** Brownian motion causes random thermal fluctuations that bring particles into contact. Smaller particles diffuse faster, leading to the inverse dependence on particle size.

**Dominant regime:** Sub-micrometer particles (< 1 μm).

### 3. Differential Sedimentation

Based on gravity-driven collisions due to differences in settling velocities:

$$\beta_{ds}(v_i, v_j) = \frac{\pi g}{72 \mu} |\rho_i - \rho_w| (d_i + d_j)^3 |d_i - d_j|$$

where:
- $g = 9.81$ m/s²: Gravitational acceleration
- $\rho_i$: Effective density of particle $i$ [kg/m³]
- $\rho_w$: Water density [kg/m³]

**Physical interpretation:** Larger/denser particles settle faster than smaller/lighter ones, creating relative motion that leads to collisions. The kernel is proportional to the difference in settling velocities.

**Dominant regime:** Large particles (> 40 μm).

---

## Breakage Kernel $S(v)$

The breakage rate represents the frequency at which particles of volume $v$ break due to hydrodynamic stress:

$$S(v_i) = k_b \left(\frac{\varepsilon \, d_i^2}{\sigma_f}\right)^n \exp\left(-\frac{C \, \sigma_f}{\rho_f \, \varepsilon^{2/3} \, d_i^{2/3}}\right)$$

where:
- $k_b$: Breakage rate constant [1/s] (typically $10^{-3}$ to $10^{-2}$)
- $n$: Breakage exponent [-] (typically 0.5 to 2.0)
- $C$: Breakage constant [-] (typically 0.5 to 2.0)
- $\sigma_f$: Floc binding strength [N/m²]
- $\rho_f$: Effective floc density [kg/m³]
- $\varepsilon$: Turbulent energy dissipation rate [m²/s³]
- $d_i$: Particle diameter [m]

### Physical Interpretation

The breakage kernel balances two competing factors:

1. **Power-law term** $\left(\frac{\varepsilon \, d_i^2}{\sigma_f}\right)^n$:
   - Represents the ratio of hydrodynamic stress to floc strength
   - Larger flocs experience stronger turbulent stresses
   - Breakage increases with turbulent intensity $\varepsilon$

2. **Exponential term** $\exp\left(-\frac{C \, \sigma_f}{\rho_f \, \varepsilon^{2/3} \, d_i^{2/3}}\right)$:
   - Represents the resistance to breakage
   - Stronger flocs (higher $\sigma_f$) are more resistant
   - Very small flocs are highly resistant due to cohesive forces

**Result:** An optimal floc size exists where aggregation and breakage balance, leading to steady-state distributions.

---

## Daughter Particle Distribution $\Gamma(v|v_{parent})$

For **binary breakage** (one particle breaks into two fragments):

$$\Gamma(v|v_{parent}) = \begin{cases}
\frac{2}{v_{parent}} & \text{if } v \leq v_{parent} \\
0 & \text{if } v > v_{parent}
\end{cases}$$

**Normalization:**
$$\int_0^{v_{parent}} \Gamma(v|v_{parent}) \, dv = 2$$

This ensures that one parent particle produces two daughter particles.

**Equal-sized daughters assumption:**
If breakage produces two equal-sized daughters:
$$v_{daughter} = \frac{v_{parent}}{2}$$

---

## Floc Properties and Fractal Structure

### Effective Floc Density

Flocs have a fractal structure, leading to size-dependent effective density:

$$\rho_{eff}(d) = \rho_w + (\rho_{primary} - \rho_w) \left(\frac{d}{d_0}\right)^{D_f - 3}$$

where:
- $\rho_w$: Water density [kg/m³]
- $\rho_{primary}$: Density of primary particles [kg/m³]
- $d_0$: Primary particle diameter [m]
- $D_f$: Fractal dimension [-] (typically 1.8-2.5 for flocs)

**Physical interpretation:** Larger flocs are more porous and have lower effective density. For $D_f < 3$, density decreases as size increases.

### Floc Mass

$$m(d) = \rho_{eff}(d) \cdot \frac{\pi d^3}{6}$$

---

## Hydrodynamic Parameters

### Relationship between Shear Rate and Dissipation Rate

In homogeneous isotropic turbulence:

$$G = \sqrt{\frac{\varepsilon}{\nu}}$$

or equivalently:

$$\varepsilon = G^2 \nu$$

where:
- $G$: Mean shear rate [1/s]
- $\varepsilon$: Turbulent energy dissipation rate [m²/s³]
- $\nu = \mu/\rho$: Kinematic viscosity [m²/s]

### Kolmogorov Scales

**Kolmogorov length scale** (smallest eddy size):

$$\eta = \left(\frac{\nu^3}{\varepsilon}\right)^{1/4}$$

**Kolmogorov time scale**:

$$\tau_{\eta} = \left(\frac{\nu}{\varepsilon}\right)^{1/2}$$

**Kolmogorov velocity scale**:

$$u_{\eta} = (\nu \varepsilon)^{1/4}$$

---

## Discretized Form (Method of Classes)

For numerical solution, the PBE is discretized into size classes (bins):

$$\frac{dn_i}{dt} = B_{agg,i} - D_{agg,i} + B_{break,i} - D_{break,i}$$

where $n_i$ is the number density in bin $i$ with representative volume $v_i$.

### Aggregation Birth (Bin i)

$$B_{agg,i} = \frac{1}{2} \sum_{j,k: v_j + v_k \approx v_i} \beta(v_j, v_k) \, n_j \, n_k$$

### Aggregation Death (Bin i)

$$D_{agg,i} = n_i \sum_{j=1}^{N_{bins}} \beta(v_i, v_j) \, n_j$$

### Breakage Birth (Bin i)

For binary breakage:

$$B_{break,i} = 2 \sum_{j: v_j \approx 2v_i} S(v_j) \, n_j$$

### Breakage Death (Bin i)

$$D_{break,i} = S(v_i) \, n_i$$

---

## Moments of the Distribution

The $k$-th moment of the particle size distribution:

$$M_k = \int_0^{\infty} v^k \, n(v) \, dv$$

**Important moments:**
- $M_0$: Total number concentration [#/m³]
- $M_1$: Total particle volume per unit volume (volume fraction) [-]
- $M_2$: Total surface area per unit volume [m²/m³]
- $M_3$: Related to mass concentration

### Volume-Weighted Mean Diameter

$$d_{43} = \frac{\sum_i d_i^4 n_i}{\sum_i d_i^3 n_i}$$

This is the **Sauter mean diameter**, commonly used to characterize floc size distributions.

---

## Conservation Properties

### Volume Conservation

In aggregation-breakage processes (without growth or dissolution), total particle volume should be conserved:

$$\frac{d}{dt} \int_0^{\infty} v \, n(v) \, dv = 0$$

**Physically:** Mass is conserved; aggregation and breakage only redistribute the size distribution without changing total mass.

### Number Conservation (Aggregation Only)

For pure aggregation without breakage:

$$\frac{d}{dt} \int_0^{\infty} n(v) \, dv = -\frac{1}{2} \int_0^{\infty} \int_0^{\infty} \beta(v_i, v_j) \, n(v_i) \, n(v_j) \, dv_i \, dv_j < 0$$

Total number decreases as particles aggregate.

---

## Typical Parameter Values

### Water Properties (20°C)

| Parameter | Symbol | Value | Units |
|-----------|--------|-------|-------|
| Water density | $\rho_w$ | 998 | kg/m³ |
| Dynamic viscosity | $\mu$ | $1.002 \times 10^{-3}$ | Pa·s |
| Kinematic viscosity | $\nu$ | $1.004 \times 10^{-6}$ | m²/s |
| Temperature | $T$ | 293.15 | K |

### Floc Properties

| Parameter | Symbol | Typical Range | Units |
|-----------|--------|---------------|-------|
| Collision efficiency | $\alpha$ | 0.01 - 0.5 | - |
| Fractal dimension | $D_f$ | 1.8 - 2.5 | - |
| Binding strength | $\sigma_f$ | $10^2$ - $10^4$ | N/m² |
| Primary particle size | $d_0$ | 0.1 - 10 | μm |

### Flocculation Conditions

| Parameter | Symbol | Typical Range | Units | Application |
|-----------|--------|---------------|-------|-------------|
| Shear rate | $G$ | 10 - 100 | 1/s | Slow mixing |
| | | 100 - 1000 | 1/s | Rapid mixing |
| Dissipation rate | $\varepsilon$ | $10^{-4}$ - $10^{-1}$ | m²/s³ | Flocculation |
| Retention time | $t$ | 10 - 1800 | s | Process design |

---

## References

1. **Smoluchowski, M.** (1917). "Versuch einer mathematischen Theorie der Koagulationskinetik kolloider Lösungen." *Z. Phys. Chem.*, 92, 129-168.

2. **Saffman, P.G., & Turner, J.S.** (1956). "On the collision of drops in turbulent clouds." *J. Fluid Mech.*, 1(1), 16-30.

3. **Pandya, J.D., & Spielman, L.A.** (1982). "Floc breakage in agitated suspensions." *J. Colloid Interface Sci.*, 90(2), 517-531.

4. **Camp, T.R., & Stein, P.C.** (1943). "Velocity gradients and internal work in fluid motion." *J. Boston Soc. Civ. Eng.*, 30, 219-237.

5. **Li, D., & Ganczarczyk, J.** (1989). "Fractal geometry of particle aggregates generated in water and wastewater treatment processes." *Environ. Sci. Technol.*, 23(11), 1385-1389.

6. **Zhan, M., et al.** (2021). "Numerical simulation of mechanical flocculation in water treatment." *Environ. Eng. Sci.*

7. **Zhan, M., et al.** (2023). "Numerical simulation of a mechanical flocculation process with multi-chambers in series." *Water Sci. Technol.*, 87(8), 1945-1960.

---

**Document Version:** 1.0
**Last Updated:** 2025-11-06
