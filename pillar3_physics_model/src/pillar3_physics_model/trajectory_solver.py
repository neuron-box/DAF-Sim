"""
Trajectory Analysis Solver for Bubble-Particle Collision

This module implements the trajectory analysis method to determine if a particle
approaching a bubble will successfully collide with it, considering both
hydrodynamic and DLVO forces.

References:
    - Han et al. (2001). Water Science and Technology, 43(8), 139-144.
    - Han, M., & Lawler, D. F. (1991). Journal AWWA, 83(3), 71-77.
"""

import math
from typing import Tuple, Optional
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp

from .particle import Particle
from .bubble import Bubble
from .fluid_properties import FluidProperties
from .dlvo_forces import DLVOForces
from .hydrodynamic_forces import HydrodynamicForces


@dataclass
class TrajectoryResult:
    """
    Result of a trajectory analysis.

    Attributes:
        collision_occurred: Whether the particle collided with the bubble
        initial_offset: Initial radial offset of particle [m]
        final_separation: Final surface-to-surface separation [m]
        trajectory_r: Array of radial positions [m]
        trajectory_theta: Array of angular positions [rad]
        trajectory_time: Array of time values [s]
        success: Whether the integration completed successfully
        message: Status message
    """
    collision_occurred: bool
    initial_offset: float
    final_separation: float
    trajectory_r: np.ndarray
    trajectory_theta: np.ndarray
    trajectory_time: np.ndarray
    success: bool
    message: str


class TrajectorySolver:
    """
    Solve particle trajectory in the flow field around a rising bubble.

    The particle motion is governed by:
    m * dv/dt = F_hydrodynamic + F_DLVO + F_gravity

    In spherical coordinates (r, θ) with origin at bubble center:
    - r: radial distance from bubble center
    - θ: angle from vertical axis (0 = top, π = bottom)
    """

    def __init__(self, particle: Particle, bubble: Bubble, fluid: FluidProperties):
        """
        Initialize trajectory solver.

        Args:
            particle: Particle object
            bubble: Bubble object
            fluid: FluidProperties object
        """
        self.particle = particle
        self.bubble = bubble
        self.fluid = fluid

        self.dlvo = DLVOForces(particle, bubble, fluid)
        self.hydro = HydrodynamicForces(particle, bubble, fluid)

        # Effective mass including added mass
        C_m = self.hydro.added_mass_coefficient()
        self.m_eff = self.particle.mass * (1.0 + C_m * fluid.density / particle.density)

    def equations_of_motion(self, t: float, state: np.ndarray) -> np.ndarray:
        """
        Equations of motion in spherical coordinates.

        State vector: [r, θ, v_r, v_θ]
        - r: radial position [m]
        - θ: angular position [rad]
        - v_r: radial velocity [m/s]
        - v_θ: angular velocity [rad/s] (not tangential velocity)

        Returns:
            Time derivatives: [dr/dt, dθ/dt, dv_r/dt, dv_θ/dt]
        """
        r, theta, v_r, v_theta_angular = state

        # Convert angular velocity to tangential velocity
        v_theta = r * v_theta_angular

        # Calculate surface-to-surface separation
        h = r - self.particle.radius - self.bubble.radius

        # Ensure minimum separation to avoid singularities
        if h < 1e-10:
            h = 1e-10

        # === DLVO Forces ===
        # DLVO force acts radially
        F_dlvo_radial = self.dlvo.total_dlvo_force(h)

        # === Hydrodynamic Forces ===
        v_particle = (v_r, v_theta)
        F_drag_r, F_drag_theta = self.hydro.drag_force_spherical(r, theta, v_particle)

        # === Gravitational Force ===
        # Convert to spherical components
        F_g_radial = self.hydro.gravitational_force() * math.cos(theta)
        F_g_tangential = -self.hydro.gravitational_force() * math.sin(theta)

        # === Total Forces ===
        F_r_total = F_dlvo_radial + F_drag_r + F_g_radial
        F_theta_total = F_drag_theta + F_g_tangential

        # === Accelerations ===
        # In spherical coordinates:
        # a_r = dv_r/dt - r * (dθ/dt)² (centrifugal term)
        # a_θ = r * d²θ/dt² + 2 * v_r * dθ/dt (Coriolis-like term)

        # From Newton's second law:
        # m * (dv_r/dt - r * ω²) = F_r
        # m * r * dω/dt + m * 2 * v_r * ω = F_θ

        omega = v_theta_angular  # angular velocity

        # Radial acceleration
        dv_r_dt = (F_r_total / self.m_eff) + r * omega**2

        # Angular acceleration
        if abs(r) > 1e-15:  # Avoid division by zero
            domega_dt = (F_theta_total / self.m_eff - 2.0 * v_r * omega) / r
        else:
            domega_dt = 0.0

        return np.array([v_r, omega, dv_r_dt, domega_dt])

    def collision_detected(self, r: float) -> bool:
        """
        Check if particle has collided with bubble.

        Args:
            r: Current radial position [m]

        Returns:
            True if collision occurred
        """
        contact_distance = self.particle.radius + self.bubble.radius
        return r <= contact_distance * 1.01  # 1% tolerance

    def particle_escaped(self, r: float, theta: float) -> bool:
        """
        Check if particle has escaped the bubble's influence.

        Args:
            r: Current radial position [m]
            theta: Current angular position [rad]

        Returns:
            True if particle escaped
        """
        # Particle escaped if:
        # 1. It's far from the bubble (> 10 * bubble radius)
        # 2. It's moving away (θ > π/2 and increasing)
        escaped_distance = self.bubble.radius * 10.0

        return r > escaped_distance or theta > math.pi

    def solve_trajectory(self, initial_offset: float, max_time: float = 1.0) -> TrajectoryResult:
        """
        Solve particle trajectory for a given initial offset.

        The particle starts far upstream (above the bubble) with an initial
        lateral offset, and the trajectory is integrated until collision or escape.

        Args:
            initial_offset: Initial radial offset from centerline [m]
                           (perpendicular distance from bubble axis)
            max_time: Maximum integration time [s]

        Returns:
            TrajectoryResult object
        """
        # Initial position: far upstream with given offset
        # Start at r = sqrt((R_b + offset)² + (3*R_b)²)
        # and angle such that horizontal offset = initial_offset
        R_b = self.bubble.radius
        z_initial = 3.0 * R_b  # Start 3 bubble radii upstream

        r_initial = math.sqrt(initial_offset**2 + z_initial**2)
        if r_initial < R_b:
            r_initial = R_b * 1.1  # Start just outside bubble

        theta_initial = math.atan2(initial_offset, z_initial)

        # Initial velocity: particle at rest in lab frame
        # (will be accelerated by bubble's flow field)
        v_r_initial = 0.0
        v_theta_initial = 0.0

        # Initial state
        state0 = np.array([r_initial, theta_initial, v_r_initial, v_theta_initial])

        # Event functions for termination
        def collision_event(t, state):
            """Stop if collision occurs."""
            r = state[0]
            return r - (self.particle.radius + self.bubble.radius) * 1.01

        def escape_event(t, state):
            """Stop if particle escapes."""
            r, theta = state[0], state[1]
            return -1.0 if self.particle_escaped(r, theta) else 1.0

        collision_event.terminal = True
        collision_event.direction = -1  # Trigger when decreasing

        escape_event.terminal = True

        # Solve ODE
        try:
            solution = solve_ivp(
                self.equations_of_motion,
                (0, max_time),
                state0,
                method='RK45',
                events=[collision_event, escape_event],
                dense_output=True,
                max_step=1e-4,
                rtol=1e-6,
                atol=1e-9
            )

            # Extract trajectory
            r_traj = solution.y[0, :]
            theta_traj = solution.y[1, :]
            t_traj = solution.t

            # Determine if collision occurred
            final_r = r_traj[-1]
            final_separation = final_r - (self.particle.radius + self.bubble.radius)

            collision_occurred = self.collision_detected(final_r)

            message = "Collision occurred" if collision_occurred else "Particle escaped"

            return TrajectoryResult(
                collision_occurred=collision_occurred,
                initial_offset=initial_offset,
                final_separation=final_separation,
                trajectory_r=r_traj,
                trajectory_theta=theta_traj,
                trajectory_time=t_traj,
                success=solution.success,
                message=message
            )

        except Exception as e:
            # Integration failed
            return TrajectoryResult(
                collision_occurred=False,
                initial_offset=initial_offset,
                final_separation=float('inf'),
                trajectory_r=np.array([r_initial]),
                trajectory_theta=np.array([theta_initial]),
                trajectory_time=np.array([0.0]),
                success=False,
                message=f"Integration failed: {str(e)}"
            )

    def find_critical_offset(self, tolerance: float = 1e-8) -> float:
        """
        Find the critical initial offset where collision just barely occurs.

        Uses bisection method to find the maximum offset that still results in collision.

        Args:
            tolerance: Tolerance for bisection [m]

        Returns:
            Critical offset [m]
        """
        # Start with bracket: [0, 2*R_bubble]
        offset_min = 0.0
        offset_max = 2.0 * self.bubble.radius

        # First check if collision occurs at all
        result_min = self.solve_trajectory(offset_min)
        if not result_min.collision_occurred:
            return 0.0  # No collision even at centerline

        result_max = self.solve_trajectory(offset_max)
        if result_max.collision_occurred:
            # Expand upper bound
            while result_max.collision_occurred:
                offset_max *= 2.0
                result_max = self.solve_trajectory(offset_max)

        # Bisection
        while (offset_max - offset_min) > tolerance:
            offset_mid = (offset_min + offset_max) / 2.0
            result_mid = self.solve_trajectory(offset_mid)

            if result_mid.collision_occurred:
                offset_min = offset_mid
            else:
                offset_max = offset_mid

        return (offset_min + offset_max) / 2.0
