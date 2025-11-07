# PBM Equation Shell v1.1.0 - Geometric Grid Support

## Summary of Changes

This update addresses critical mass conservation issues identified in PR code review and implements proper support for geometric particle volume grids.

## Critical Fixes

### 1. **CRITICAL: Fixed Mass Conservation for Geometric Grids**
- **Issue**: Original implementation assumed linear volume grids (v_i = i·Δv), but examples used geometric grids (v_i = v_0·r^i)
- **Impact**: Mass was NOT conserved during aggregation, violating fundamental physics
- **Fix**: Implemented fixed pivot technique (Kumar & Ramkrishna, 1996) to properly distribute aggregates between adjacent bins
- **Reference**: PR Comment #1 (Critical severity)

### 2. **HIGH: Fixed False Positive in Mass Conservation Test**
- **Issue**: Test used unrealistically small beta values (1e-18), hiding the conservation error
- **Impact**: Critical bugs were not detected by test suite
- **Fix**: Updated tests to use realistic beta values (1e-12) and proper grid configurations
- **Reference**: PR Comment #2 (High severity)

## API Changes

### Breaking Changes

#### `PopulationBalanceModel.__init__()`
**Old signature:**
```python
PopulationBalanceModel(n_classes: int, validate_inputs: bool = True)
```

**New signature:**
```python
PopulationBalanceModel(
    n_classes: int,
    particle_volumes: np.ndarray,
    validate_inputs: bool = True,
    grid_type: Optional[Literal['linear', 'geometric']] = None
)
```

**Migration:**
```python
# Old code (BROKEN):
pbm = PopulationBalanceModel(n_classes=5)

# New code (CORRECT):
particle_volumes = np.array([1e-18, 2e-18, 4e-18, 8e-18, 16e-18])
pbm = PopulationBalanceModel(n_classes=5, particle_volumes=particle_volumes)
```

#### `validate_mass_conservation()`
**Old signature:**
```python
validate_mass_conservation(N, dndt, particle_volumes, tolerance=1e-10)
```

**New signature:**
```python
validate_mass_conservation(N, dndt, tolerance=1e-10)
```

Particle volumes are now stored in the class, so no need to pass them again.

## New Features

### 1. Automatic Grid Type Detection
The class automatically detects whether your grid is linear or geometric:

```python
# Geometric grid (auto-detected)
vols_geo = np.array([1e-18, 2e-18, 4e-18, 8e-18, 16e-18])
pbm_geo = PopulationBalanceModel(5, vols_geo)
print(pbm_geo.grid_type)  # 'geometric'
print(pbm_geo.grid_ratio)  # 2.0

# Linear grid (auto-detected)
vols_lin = np.array([1e-18, 2e-18, 3e-18, 4e-18, 5e-18])
pbm_lin = PopulationBalanceModel(5, vols_lin)
print(pbm_lin.grid_type)  # 'linear'
```

### 2. Fixed Pivot Technique for Geometric Grids
When two particles aggregate on a geometric grid, the resulting volume typically doesn't match any grid point exactly. The fixed pivot technique distributes the mass between two adjacent bins to conserve both particle number and volume:

```
v_new = v_j + v_k
If v_i ≤ v_new < v_{i+1}:
  ξ = (v_{i+1} - v_new) / (v_{i+1} - v_i)  → fraction to bin i
  η = (v_new - v_i) / (v_{i+1} - v_i)      → fraction to bin i+1
```

This ensures exact mass conservation for particles within the grid.

### 3. Explicit Grid Type Specification
You can force a specific grid type interpretation:

```python
pbm = PopulationBalanceModel(
    n_classes=5,
    particle_volumes=vols,
    grid_type='geometric'  # or 'linear'
)
```

## Important Notes on Mass Conservation

### Expected Behavior
- **Geometric grids with fixed pivot**: Mass is conserved EXACTLY for particles that stay within the grid
- **Finite grid overflow**: Particles that aggregate beyond v_max are placed in the largest bin, causing mass loss
- **Linear grids**: Poor mass conservation because v_i + v_j rarely equals any v_k exactly

### Recommendations
1. Use **geometric grids** whenever possible (better mass conservation)
2. Choose grid range to minimize overflow:
   - Max volume should be >> largest expected aggregate
   - For conservative design: v_max ≥ 2 * v_{n-2} (handles largest collisions)
3. Monitor mass conservation in your simulations
4. Use grid ratio r ≈ 2 for best trade-off between coverage and resolution

### Example of Good Grid Design
```python
# Design for particles up to 1000 μm
n_classes = 15
d_min, d_max = 1e-6, 2000e-6  # 1 μm to 2mm (safety factor of 2)
particle_diameters = np.logspace(np.log10(d_min), np.log10(d_max), n_classes)
particle_volumes = (4/3) * np.pi * (particle_diameters/2)**3

pbm = PopulationBalanceModel(n_classes, particle_volumes)
```

## Documentation of Limitations (PR Comment #3 - Medium)

Added clear documentation of `construct_beta_matrix()` limitations:
- Assumes constant shear rate (no spatial/temporal variation)
- Does not include hydrodynamic interactions
- Uses point-particle approximation (no fractal dimensions)
- Does not account for Hamaker constant

These are acceptable for v1.0 as the function is a helper - users can provide their own beta matrices.

## Test Suite Improvements

- Added 10+ new tests specifically for geometric grids
- Tests now use realistic beta values (exposed the mass conservation bug)
- Added tests for grid detection, distribution factors, and overflow handling
- Mass conservation tests now properly account for finite grid limitations
- All 31 tests pass

## Backward Compatibility

**NOT backward compatible** due to required `particle_volumes` parameter. This is necessary to fix the critical mass conservation bug.

## Upgrade Guide

1. Add `particle_volumes` array to all `PopulationBalanceModel` instantiations
2. Remove `particle_volumes` argument from `validate_mass_conservation()` calls
3. Verify your grid type (geometric recommended)
4. Check mass conservation in your workflows

## References

- Kumar, S., & Ramkrishna, D. (1996). "On the solution of population balance equations by discretization—I. A fixed pivot technique." Chemical Engineering Science, 51(8), 1311-1332.
- Abreu, D. A. d. S., et al. (2021). "Integration of first principle models and machine learning in a modeling framework: An application to flocculation." Chemical Engineering Science.

## Version

- Previous: v1.0.0 (BROKEN mass conservation)
- Current: v1.1.0 (FIXED with geometric grid support)
