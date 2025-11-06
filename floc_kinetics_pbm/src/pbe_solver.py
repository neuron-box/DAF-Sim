"""
Population Balance Equation (PBE) solver for floc kinetics.

This module implements the numerical solution of the PBE for floc aggregation
and breakage in a spatially discretized domain (e.g., CFD cells).
"""

import numpy as np
from typing import Optional, Dict, Tuple, Callable
from scipy.integrate import odeint
from .kernels import FlocKernels
from .properties import FlocProperties


class PopulationBalanceModel:
    """
    Solves the Population Balance Equation for floc size distribution.

    The PBE describes the evolution of particle number density n(v,x,t) where:
    - v: particle volume
    - x: spatial position
    - t: time

    Discretized form (method of classes):
    dn_i/dt = B_agg,i - D_agg,i + B_break,i - D_break,i + transport_i

    where i represents the size class.
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

        # Bin widths
        self.dv = np.zeros(n_bins)
        self.dv[0] = self.volumes[1] - self.volumes[0]
        for i in range(1, n_bins - 1):
            self.dv[i] = 0.5 * (self.volumes[i + 1] - self.volumes[i - 1])
        self.dv[-1] = self.volumes[-1] - self.volumes[-2]

        # Precompute aggregation kernel matrix for efficiency
        self.beta_matrix = None
        self.breakage_rates = None

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
        Calculate aggregation birth term for bin i.

        B_agg,i = (1/2) * Σ_{j,k: v_j+v_k=v_i} β(j,k) * n_j * n_k

        Args:
            n: Number density array [#/m³]
            i: Bin index

        Returns:
            Aggregation birth rate [#/m³/s]
        """
        if self.beta_matrix is None:
            return 0.0

        birth = 0.0
        v_target = self.volumes[i]

        # Sum over all pairs that aggregate to form particles in bin i
        for j in range(self.n_bins):
            for k in range(j, self.n_bins):  # k >= j to avoid double counting
                v_sum = self.volumes[j] + self.volumes[k]

                # Check if aggregation produces particle in bin i
                if abs(v_sum - v_target) < 0.5 * self.dv[i]:
                    if j == k:
                        birth += 0.5 * self.beta_matrix[j, k] * n[j] * n[k]
                    else:
                        birth += self.beta_matrix[j, k] * n[j] * n[k]

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

        # Precompute kernels for efficiency
        if recompute_kernels or self.beta_matrix is None:
            self.precompute_kernels(G)

        # Solve ODE
        n_solution = odeint(self.pbe_rhs, n_initial, t_span, args=(G,))

        return t_span, n_solution

    def moments(self, n: np.ndarray, k: int = 0) -> float:
        """
        Calculate k-th moment of the distribution.

        M_k = Σ_i v_i^k * n_i * Δv_i

        Args:
            n: Number density array [#/m³]
            k: Moment order

        Returns:
            k-th moment
        """
        return np.sum(self.volumes**k * n * self.dv)

    def mean_diameter(self, n: np.ndarray) -> float:
        """
        Calculate volume-weighted mean diameter (d_43).

        d_43 = M_4 / M_3 where M_k is the k-th moment based on diameter.

        Args:
            n: Number density array [#/m³]

        Returns:
            Mean diameter [m]
        """
        M3 = np.sum(self.diameters**3 * n * self.dv)
        M4 = np.sum(self.diameters**4 * n * self.dv)

        if M3 > 0:
            return M4 / M3
        else:
            return 0.0

    def total_number_concentration(self, n: np.ndarray) -> float:
        """
        Calculate total number concentration.

        N_total = Σ_i n_i * Δv_i

        Args:
            n: Number density array [#/m³]

        Returns:
            Total number concentration [#/m³]
        """
        return self.moments(n, k=0)

    def total_volume_fraction(self, n: np.ndarray) -> float:
        """
        Calculate total volume fraction of particles.

        φ = Σ_i v_i * n_i * Δv_i

        Args:
            n: Number density array [#/m³]

        Returns:
            Volume fraction [-]
        """
        return self.moments(n, k=1)

    def summary_statistics(self, n: np.ndarray) -> Dict[str, float]:
        """
        Calculate summary statistics of the particle size distribution.

        Args:
            n: Number density array [#/m³]

        Returns:
            Dictionary with statistics:
                - total_number: Total number concentration [#/m³]
                - volume_fraction: Total volume fraction [-]
                - mean_diameter: Volume-weighted mean diameter [m]
                - d10: 10th percentile diameter [m]
                - d50: Median diameter [m]
                - d90: 90th percentile diameter [m]
        """
        # Cumulative distribution
        cumulative = np.cumsum(n * self.dv)
        if cumulative[-1] > 0:
            cumulative = cumulative / cumulative[-1]
        else:
            cumulative = np.zeros_like(cumulative)

        # Percentiles
        d10 = np.interp(0.1, cumulative, self.diameters)
        d50 = np.interp(0.5, cumulative, self.diameters)
        d90 = np.interp(0.9, cumulative, self.diameters)

        return {
            'total_number': self.total_number_concentration(n),
            'volume_fraction': self.total_volume_fraction(n),
            'mean_diameter': self.mean_diameter(n),
            'd10': d10,
            'd50': d50,
            'd90': d90
        }
