# Geometric Grid Implementation Notes

## Problem Statement

For geometric grids where particle volumes follow: v_i = v_0 × r^i

When two particles collide (j + k), the resulting volume v_new = v_j + v_k typically does NOT match any grid point exactly.

## Solution: Fixed Pivot Technique (Kumar & Ramkrishna 1996)

### Key Concept
Distribute the newly formed particle between two adjacent bins to conserve mass.

### Algorithm

1. **Find target bin:** For v_new = v_j + v_k, find index i where:
   ```
   v_i ≤ v_new < v_{i+1}
   ```

2. **Calculate distribution factors:**

   To conserve both number and volume, use:
   ```
   ξ = (v_{i+1} - v_new) / (v_{i+1} - v_i)
   η = (v_new - v_i) / (v_{i+1} - v_i)
   ```
   where ξ + η = 1

3. **Distribute births:**
   - Fraction ξ goes to bin i
   - Fraction η goes to bin i+1

### Special Cases

1. **v_new ≥ v_max:** All mass to largest bin (loss term)
2. **v_new < v_min:** Shouldn't happen in practice
3. **Linear grid:** ξ calculation simplifies, reduces to original implementation

### Modified PBE for Geometric Grids

Birth by aggregation becomes:
```
Birth_i = (1/2) × Σ_{j,k: v_i ≤ v_j+v_k < v_{i+1}} ξ_{jk,i} × β_{jk} × α_{jk} × N_j × N_k
        + (1/2) × Σ_{j,k: v_{i-1} ≤ v_j+v_k < v_i} η_{jk,i} × β_{jk} × α_{jk} × N_j × N_k
```

where:
- First term: particles landing in [v_i, v_{i+1}) with fraction ξ to bin i
- Second term: particles landing in [v_{i-1}, v_i) with fraction η to bin i

## Implementation Strategy

### Phase 1: Add Volume Grid Support
- Add `particle_volumes` parameter to `__init__`
- Auto-detect grid type (linear vs geometric)
- Store grid ratio `r` for geometric grids

### Phase 2: Implement Distribution Logic
- Helper function to find target bin for v_new
- Helper function to calculate ξ and η
- Modify aggregation birth term to distribute mass

### Phase 3: Validation
- Test mass conservation with realistic beta values
- Compare linear vs geometric grid results
- Benchmark against literature examples

## References

1. Kumar, S., & Ramkrishna, D. (1996). "On the solution of population balance equations by discretization—I. A fixed pivot technique." Chemical Engineering Science, 51(8), 1311-1332.

2. Hounslow, M. J., Ryall, R. L., & Marshall, V. R. (1988). "A discretized population balance for nucleation, growth, and aggregation." AIChE Journal, 34(11), 1821-1832.

3. Abreu, D. A. d. S., et al. (2021). "Integration of first principle models and machine learning in a modeling framework: An application to flocculation." Chemical Engineering Science.
