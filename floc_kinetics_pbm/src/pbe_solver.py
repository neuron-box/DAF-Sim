"""
Population Balance Equation (PBE) solver for floc kinetics.

This module implements the numerical solution of the PBE for floc aggregation
and breakage in a spatially discretized domain (e.g., CFD cells).

Discretization Convention:
--------------------------
The continuous PBE is discretized using the "method of classes":
- Size domain divided into n_bins
- n[i] = number density in bin i [particles/m³ of fluid]
- v[i] = representative volume of bin i [m³]
- dv[i] = width of bin i in volume space [m³]

The number density n[i] represents particles per unit volume of FLUID
(not per bin width). Therefore:
  - Total particles = Σ n[i] * dv[i] [dimensionless, particles per m³ fluid]
  - Total volume = Σ n[i] * v[i] * dv[i] [m³ particles / m³ fluid]

Note: The factor dv[i] appears because we're integrating over a discrete
representation of a continuous distribution. Think of n[i]*dv[i] as the
integral ∫_{bin i} n(v) dv.
"""

import numpy as np
from typing import Optional, Dict, Tuple, Callable, List
from scipy.integrate import odeint
from .kernels import FlocKernels
from .properties import FlocProperties


class PopulationBalanceModel:
    """
    Solves the Population Balance Equation for floc size distribution.

    The PBE describes the evolution of particle number density n(v,x,t) where:
    - v: particle volume [m³]
    - x: spatial position [m]
    - t: time [s]

    Discretized form (method of classes):
    dn_i/dt = B_agg,i - D_agg,i + B_break,i - D_break,i + transport_i

    where i represents the size class (bin).

    Conservation Property:
    ---------------------
    During aggregation and breakage (without transport), total particle
    volume should be conserved:
        d/dt[Σ v_i * n_i * dv_i] ≈ 0

    Statistics Definitions:
    ----------------------
    - mean_diameter: Sauter mean diameter d₄₃ = Σ(n·dv·d⁴)/Σ(n·dv·d³)
    - d10_number, d50_number, d90_number: Number-based percentiles
    - d10_volume, d50_volume, d90_volume: Volume-based percentiles
    - total_number: Total number concentration Σ n·dv [#/m³]
    - volume_fraction: Volume occupied by particles Σ n·v·dv [-]
    """

    def __init__(
        self,
        n_bins: int = 20,
        d_min: float = 1e-6,
        d_max: float = 1e-3,
        properties: Optional[FlocProperties] = None,
        kernels: Optional[FlocKernels] = None
    ):
        """
        Initialize the PBM solver.

        Args:
            n_bins: Number of size classes (bins)
            d_min: Minimum particle diameter [m]
            d_max: Maximum particle diameter [m]
            properties: FlocProperties object
            kernels: FlocKernels object
        """
        self.n_bins = n_bins
        self.d_min = d_min
        self.d_max = d_max

        # Initialize properties and kernels
        self.props = properties if properties is not None else FlocProperties()
        self.kernels = kernels if kernels is not None else FlocKernels(self.props)

        # Create size bins (logarithmic spacing is typical for PBM)
        self.diameters = np.logspace(np.log10(d_min), np.log10(d_max), n_bins)
        self.volumes = self.kernels.volume_from_diameter(self.diameters)

        # Bin widths in volume space
        self.dv = np.zeros(n_bins)
        self.dv[0] = self.volumes[1] - self.volumes[0]
        for i in range(1, n_bins - 1):
            self.dv[i] = 0.5 * (self.volumes[i + 1] - self.volumes[i - 1])
        self.dv[-1] = self.volumes[-1] - self.volumes[-2]

        # Precompute aggregation kernel matrix for efficiency
        self.beta_matrix = None
        self.breakage_rates = None

        # Aggregation mapping for O(n²) performance
        self.agg_map = None

    def precompute_kernels(self, G: float):
        """
        Precompute aggregation kernel matrix and breakage rates.

        This significantly speeds up the time integration by avoiding
        repeated kernel calculations.

        Args:
            G: Shear rate [1/s] (assumed constant for this computation)
        """
        n = self.n_bins

        # Aggregation kernel matrix β(i,j)
        self.beta_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                self.beta_matrix[i, j] = self.kernels.beta_total(
                    self.diameters[i], self.diameters[j], G
                )

        # Breakage rates S(i)
        self.breakage_rates = np.zeros(n)
        for i in range(n):
            self.breakage_rates[i] = self.kernels.breakage_rate(self.diameters[i], G)

        # Build aggregation mapping for optimized computation
        self._build_aggregation_map()

    def _build_aggregation_map(self):
        """
        Build a mapping of (j, k) pairs to target bins for aggregation.

        This pre-computation reduces aggregation_birth from O(n²) to O(1)
        per bin, reducing overall PBE RHS from O(n³) to O(n²).

        The mapping stores, for each target bin i, a list of (j, k, weight)
        tuples where:
        - j, k are the indices of particles that aggregate
        - weight is 0.5 if j==k (to avoid double counting), 1.0 otherwise
        """
        self.agg_map = {i: [] for i in range(self.n_bins)}

        for j in range(self.n_bins):
            for k in range(j, self.n_bins):
                v_sum = self.volumes[j] + self.volumes[k]

                # Find target bin - use closest bin
                if v_sum <= self.volumes[0]:
                    target = 0
                elif v_sum >= self.volumes[-1]:
                    target = self.n_bins - 1
                else:
                    # Binary search for closest bin
                    idx = np.searchsorted(self.volumes, v_sum)
                    if idx == 0:
                        target = 0
                    elif idx >= self.n_bins:
                        target = self.n_bins - 1
                    else:
                        # Check which bin is closer
                        if abs(v_sum - self.volumes[idx-1]) < abs(v_sum - self.volumes[idx]):
                            target = idx - 1
                        else:
                            target = idx

                # Add to mapping
                weight = 0.5 if j == k else 1.0
                self.agg_map[target].append((j, k, weight))

    def find_bin_index(self, volume: float) -> int:
        """
        Find the bin index for a given volume.

        Args:
            volume: Particle volume [m³]

        Returns:
            Bin index (0 to n_bins-1)
        """
        if volume <= self.volumes[0]:
            return 0
        if volume >= self.volumes[-1]:
            return self.n_bins - 1

        # Binary search
        idx = np.searchsorted(self.volumes, volume)
        return min(idx, self.n_bins - 1)

    def aggregation_birth(self, n: np.ndarray, i: int) -> float:
        """
        Calculate aggregation birth term for bin i (OPTIMIZED VERSION).

        B_agg,i = (1/2) * Σ_{j,k: v_j+v_k≈v_i} β(j,k) * n_j * n_k

        Uses pre-computed aggregation mapping for O(1) complexity per bin.

        Args:
            n: Number density array [#/m³]
            i: Bin index

        Returns:
            Aggregation birth rate [#/m³/s]
        """
        if self.beta_matrix is None or self.agg_map is None:
            return 0.0

        birth = 0.0
        for j, k, weight in self.agg_map[i]:
            birth += weight * self.beta_matrix[j, k] * n[j] * n[k]

        return birth

    def aggregation_death(self, n: np.ndarray, i: int) -> float:
        """
        Calculate aggregation death term for bin i.

        D_agg,i = n_i * Σ_j β(i,j) * n_j

        Args:
            n: Number density array [#/m³]
            i: Bin index

        Returns:
            Aggregation death rate [#/m³/s]
        """
        if self.beta_matrix is None:
            return 0.0

        death = n[i] * np.sum(self.beta_matrix[i, :] * n)
        return death

    def breakage_birth(self, n: np.ndarray, i: int) -> float:
        """
        Calculate breakage birth term for bin i.

        B_break,i = Σ_{j: v_j>v_i} S(j) * Γ(v_i|v_j) * n_j

        Assumes binary breakage into two equal daughters.

        Args:
            n: Number density array [#/m³]
            i: Bin index

        Returns:
            Breakage birth rate [#/m³/s]
        """
        if self.breakage_rates is None:
            return 0.0

        birth = 0.0
        v_i = self.volumes[i]

        # Particles that break to form particles in bin i
        # Binary breakage: parent volume = 2 * v_i
        v_parent = 2.0 * v_i
        j_parent = self.find_bin_index(v_parent)

        if j_parent < self.n_bins:
            # Binary breakage produces 2 daughters
            birth = 2.0 * self.breakage_rates[j_parent] * n[j_parent]

        return birth

    def breakage_death(self, n: np.ndarray, i: int) -> float:
        """
        Calculate breakage death term for bin i.

        D_break,i = S(i) * n_i

        Args:
            n: Number density array [#/m³]
            i: Bin index

        Returns:
            Breakage death rate [#/m³/s]
        """
        if self.breakage_rates is None:
            return 0.0

        death = self.breakage_rates[i] * n[i]
        return death

    def pbe_rhs(self, n: np.ndarray, t: float, G: float) -> np.ndarray:
        """
        Right-hand side of the PBE for time integration.

        dn/dt = B_agg - D_agg + B_break - D_break

        Args:
            n: Number density array [#/m³]
            t: Time [s] (not used but required by odeint)
            G: Shear rate [1/s]

        Returns:
            Time derivative of number density [#/m³/s]
        """
        # Validate input
        if np.any(n < 0):
            # Allow small negative values due to numerical errors
            n = np.maximum(n, 0)

        dndt = np.zeros(self.n_bins)

        for i in range(self.n_bins):
            # Aggregation terms
            B_agg = self.aggregation_birth(n, i)
            D_agg = self.aggregation_death(n, i)

            # Breakage terms
            B_break = self.breakage_birth(n, i)
            D_break = self.breakage_death(n, i)

            # Total rate of change
            dndt[i] = B_agg - D_agg + B_break - D_break

        return dndt

    def solve(
        self,
        n_initial: np.ndarray,
        t_span: np.ndarray,
        G: float,
        recompute_kernels: bool = True
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Solve the PBE for given initial conditions.

        Args:
            n_initial: Initial number density for each bin [#/m³]
            t_span: Time points for solution [s]
            G: Shear rate [1/s] (assumed constant)
            recompute_kernels: Whether to recompute kernel matrix

        Returns:
            Tuple of (t, n) where:
                - t: Time array [s]
                - n: Number density array [#/m³] (n_times × n_bins)
        """
        if len(n_initial) != self.n_bins:
            raise ValueError(f"Initial condition must have {self.n_bins} bins")

        if np.any(n_initial < 0):
            raise ValueError("Initial number density cannot be negative")

        # Precompute kernels for efficiency
        if recompute_kernels or self.beta_matrix is None:
            self.precompute_kernels(G)

        # Solve ODE
        n_solution = odeint(self.pbe_rhs, n_initial, t_span, args=(G,))

        return t_span, n_solution

    def moments(self, n: np.ndarray, k: int = 0) -> float:
        """
        Calculate k-th moment of the distribution in volume space.

        M_k = Σ_i v_i^k * n_i * dv_i

        Args:
            n: Number density array [#/m³]
            k: Moment order

        Returns:
            k-th moment
        """
        return np.sum(self.volumes**k * n * self.dv)

    def mean_diameter(self, n: np.ndarray) -> float:
        """
        Calculate Sauter mean diameter (d₄₃) - volume-weighted mean.

        CORRECTED FORMULA:
        d₄₃ = Σ(n_i * dv_i * d_i⁴) / Σ(n_i * dv_i * d_i³)

        This is the volume-weighted mean diameter, commonly used in
        particle characterization and flocculation studies.

        Args:
            n: Number density array [#/m³]

        Returns:
            Sauter mean diameter [m]
        """
        # Weight by particle count in each bin: n[i] * dv[i]
        particle_count = n * self.dv

        M3 = np.sum(particle_count * self.diameters**3)
        M4 = np.sum(particle_count * self.diameters**4)

        if M3 > 0:
            return M4 / M3
        else:
            return 0.0

    def total_number_concentration(self, n: np.ndarray) -> float:
        """
        Calculate total number concentration.

        CORRECTED FORMULA:
        N_total = Σ_i n_i * dv_i [#/m³]

        This integrates the number density over all size classes.

        Args:
            n: Number density array [#/m³]

        Returns:
            Total number concentration [#/m³]
        """
        return np.sum(n * self.dv)

    def total_volume_fraction(self, n: np.ndarray) -> float:
        """
        Calculate total volume fraction of particles.

        φ = Σ_i v_i * n_i * dv_i [m³ particles / m³ fluid]

        Args:
            n: Number density array [#/m³]

        Returns:
            Volume fraction [-]
        """
        return self.moments(n, k=1)

    def percentiles_number_based(self, n: np.ndarray) -> Tuple[float, float, float]:
        """
        Calculate NUMBER-BASED percentiles.

        These represent particle diameters where X% of particles are smaller.
        E.g., d50_number means 50% of particles (by count) are smaller.

        Args:
            n: Number density array [#/m³]

        Returns:
            Tuple of (d10, d50, d90) in meters
        """
        # Particle count in each bin
        particle_count = n * self.dv

        # Cumulative number fraction
        cumulative = np.cumsum(particle_count)
        if cumulative[-1] > 0:
            cumulative = cumulative / cumulative[-1]
        else:
            cumulative = np.zeros_like(cumulative)

        # Calculate percentiles
        d10 = np.interp(0.1, cumulative, self.diameters)
        d50 = np.interp(0.5, cumulative, self.diameters)
        d90 = np.interp(0.9, cumulative, self.diameters)

        return d10, d50, d90

    def percentiles_volume_based(self, n: np.ndarray) -> Tuple[float, float, float]:
        """
        Calculate VOLUME-BASED percentiles (more common in engineering).

        These represent particle diameters where X% of total particle
        volume is in smaller particles.
        E.g., d50_volume means 50% of total volume is in particles smaller.

        Args:
            n: Number density array [#/m³]

        Returns:
            Tuple of (d10, d50, d90) in meters
        """
        # Volume of particles in each bin
        volume_in_bin = n * self.dv * self.volumes

        # Cumulative volume fraction
        cumulative = np.cumsum(volume_in_bin)
        if cumulative[-1] > 0:
            cumulative = cumulative / cumulative[-1]
        else:
            cumulative = np.zeros_like(cumulative)

        # Calculate percentiles
        d10 = np.interp(0.1, cumulative, self.diameters)
        d50 = np.interp(0.5, cumulative, self.diameters)
        d90 = np.interp(0.9, cumulative, self.diameters)

        return d10, d50, d90

    def summary_statistics(self, n: np.ndarray) -> Dict[str, float]:
        """
        Calculate summary statistics of the particle size distribution.

        CORRECTED VERSION with both number-based and volume-based percentiles.

        Args:
            n: Number density array [#/m³]

        Returns:
            Dictionary with statistics:
                - total_number: Total number concentration [#/m³]
                - volume_fraction: Total volume fraction [-]
                - mean_diameter: Sauter mean diameter d₄₃ [m]
                - d10_number, d50_number, d90_number: Number-based percentiles [m]
                - d10_volume, d50_volume, d90_volume: Volume-based percentiles [m]
        """
        # Get number-based percentiles
        d10_num, d50_num, d90_num = self.percentiles_number_based(n)

        # Get volume-based percentiles
        d10_vol, d50_vol, d90_vol = self.percentiles_volume_based(n)

        return {
            'total_number': self.total_number_concentration(n),
            'volume_fraction': self.total_volume_fraction(n),
            'mean_diameter': self.mean_diameter(n),
            'd10_number': d10_num,
            'd50_number': d50_num,
            'd90_number': d90_num,
            'd10_volume': d10_vol,
            'd50_volume': d50_vol,
            'd90_volume': d90_vol,
            # For backward compatibility, use volume-based as default
            'd10': d10_vol,
            'd50': d50_vol,
            'd90': d90_vol,
        }
