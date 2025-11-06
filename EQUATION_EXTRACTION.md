# Eulerian Phase Momentum Equation - VTT (2023)

## Source Document
**Title:** Multiphase Flow Simulation of ITTC Standard Cavitator for Underwater Radiated Noise Prediction
**Authors:** Hynninen A., Viitanen V., Tanttari J., Klose R., Testa C., Martio J.
**Publication:** Journal of Marine Science and Engineering, 2023, Volume 11, Issue 4, Article 820
**DOI:** 10.3390/jmse11040820
**Institution:** VTT Technical Research Centre of Finland

## Complete Phase Momentum Equation

The Eulerian phase momentum equation as extracted from the VTT report [19, 20]:

```latex
\frac{\partial\alpha_{\phi}\rho_{\phi}U_{\phi}}{\partial t} + \nabla \cdot (\alpha_{\phi}\rho_{\phi}U_{\phi}U_{\phi}) − S_{ce,\phi}U_{\phi} + R_{\phi} + \alpha_{\phi}\nabla p − \rho_{\phi}g = M_{\phi}
```

or in display form:

$$\frac{\partial\alpha_{\phi}\rho_{\phi}U_{\phi}}{\partial t} + \nabla \cdot (\alpha_{\phi}\rho_{\phi}U_{\phi}U_{\phi}) − S_{ce,\phi}U_{\phi} + R_{\phi} + \alpha_{\phi}\nabla p − \rho_{\phi}g = M_{\phi}$$

## Term Definitions

| Mathematical Term | Physical Meaning |
|------------------|------------------|
| $\frac{\partial\alpha_{\phi}\rho_{\phi}U_{\phi}}{\partial t}$ | **Transient Term** - Rate of change of momentum with time for phase φ. Represents unsteady accumulation or depletion of momentum. |
| $\nabla \cdot (\alpha_{\phi}\rho_{\phi}U_{\phi}U_{\phi})$ | **Convection/Advection Term** - Momentum flux due to bulk fluid motion (transport of momentum by the flow itself). Represents momentum transfer due to fluid convection. |
| $S_{ce,\phi}U_{\phi}$ | **Compression-Expansion Source Term** - Momentum source/sink due to phase change, mass transfer, or compressibility effects. The subscript "ce" indicates compression-expansion phenomena. |
| $R_{\phi}$ | **Reynolds/Turbulent Stress Term** - Momentum transfer due to turbulent fluctuations and viscous stresses. Represents the divergence of the stress tensor: $R_{\phi} = \nabla \cdot \tau_{\phi}$ where $\tau_{\phi}$ includes both viscous and turbulent (Reynolds) stresses. |
| $\alpha_{\phi}\nabla p$ | **Pressure Gradient Term** - Force due to pressure variation in the flow field, weighted by the volume fraction of phase φ. Drives flow from high to low pressure regions. |
| $\rho_{\phi}g$ | **Gravitational Body Force** - Gravity force acting on the phase (note: should typically be $\alpha_{\phi}\rho_{\phi}g$ with volume fraction, but notation may vary). |
| $M_{\phi}$ | **Interfacial Momentum Transfer Term** - Net momentum exchange between phases due to interfacial forces including: drag force, lift force, virtual mass force, wall lubrication force, turbulent dispersion force, and other interfacial effects. |

## Phase Notation

- $\phi$ : Phase indicator (e.g., liquid, gas, solid particles)
- $\alpha_{\phi}$ : Volume fraction of phase φ (dimensionless, 0 ≤ α_φ ≤ 1)
- $\rho_{\phi}$ : Density of phase φ [kg/m³]
- $U_{\phi}$ : Velocity vector of phase φ [m/s]
- $p$ : Pressure field [Pa]
- $g$ : Gravitational acceleration vector [m/s²]

## Equation Form and Sign Conventions

This equation is written in a form where:
- All terms are moved to the left-hand side except the interfacial momentum transfer
- $M_{\phi}$ is isolated on the right-hand side
- Negative signs indicate forces opposing the flow direction
- The sum of interfacial forces over all phases must equal zero (action-reaction): $\sum_{\phi} M_{\phi} = 0$

## Application to DAF (Dissolved Air Flotation)

In the context of DAF simulation, this equation governs:
1. **Liquid phase (water):** Continuous phase motion with suspended particles
2. **Gas phase (microbubbles):** Rising bubble motion due to buoyancy and drag
3. **Solid phase (flocs/particles):** Particle motion and attachment to bubbles

The interfacial term $M_{\phi}$ is critical for DAF as it captures:
- Bubble-water drag interactions
- Particle-water drag
- Bubble-particle collision and attachment forces
- Three-phase interaction dynamics
