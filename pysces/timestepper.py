"""A module to easily set up and manage a simulation"""
import numpy as np
from .vortex import Vortices

__all__ = ['ExplicitEuler', 'RungeKutta2', 'RungeKutta4']

class Timestepper(object):
    """Base class for timesteppers for unsteady boundary element simulation"""

    def __init__(self, dt, Uinfty=(1,0), body=None, body_cls=None, wake=None):
        """Initialize a simulation"""
        self._dt = dt
        self._Uinfty = np.array(Uinfty)
        self._body = body
        if body is None:
            self._has_body = False
            self._bound = None
        else:
            self._has_body = True
            self._bound = body_cls(body)
        self.initialize(wake)

    def initialize(self, wake=None):
        """Initialize a timestepper

        Solve for panel strengths so that surface boundary conditions are
        satisfied, and shed a particle into the wake so that overall circulation
        is zero.

        """
        self._time = 0
        if wake is None:
            self._wake = Vortices()
        else:
            self._wake = Vortices(wake.positions, wake.strengths)
        if self._has_body:
            self._body.time = 0
            self._bound.update_strengths_unsteady(self._dt, self._Uinfty)
            self._wake.append(*self._bound.get_newly_shed())

    @property
    def time(self):
        """Current simulation time"""
        return self._time

    @property
    def bound(self):
        """Body panels used in the simulation"""
        return self._bound

    @property
    def wake(self):
        """Wake vortices used in the simulation"""
        return self._wake

    @property
    def dt(self):
        """Timestep for the simulation"""
        return self._dt

    def _wake_velocity(self, pos=None, dt=0):
        """Compute the induced velocity at each of the wake vortices

        This is the right-hand side for the timestepper that advances the
        positions of the wake vortices

        Parameters
        ----------
        pos : array, optional
            Array (shape (n,2)) of positions of wake vortices.  Default is the
            current positions of wake vortices.
        dt : float, optional
            Timestep between current simulation time, and time at which the
            velocity is to be computed (default is 0).

        Returns
        -------
        vel : array, shape (n,2)
            Induced velocities at the specified locations of wake vortices.

        Notes
        -----
        If ``pos`` is specified, the wake positions in the Timestepper object
        are updated, and the strengths of bound elements are updated as well to
        satisfy the no-flow-through boundary condition.

        """
        if pos is None:
            pos = self._wake.positions
            shed = None
        else:
            self._wake.positions = pos
            if self._has_body:
                # update body position and strengths of surface elements
                self._body.time = self._time + dt
                self._bound.update_strengths_unsteady(dt, self._Uinfty)
                shed = Vortices(*self._bound.get_newly_shed())
        vel = self._wake.induced_velocity(pos)
        vel += self._Uinfty
        if self._has_body:
            vel += self._bound.induced_velocity(pos)
            if shed:
                vel += shed.induced_velocity(pos)
        return vel

    def _update_flow(self, wake_pos, dt):
        """Update the flow with new positions of wake vortices

        Parameters
        ----------
        wake_pos : array
            The new locations of wake vortices
        dt : float
            The amount by which the time should be incremented

        Notes
        -----
        The body motion is updated to the new time, the strengths of the
        bound elements are updated to enforce the no-flow-through boundary
        condition, and a newly shed vortex is added to the wake.

        """
        self._wake.positions = wake_pos
        self._time += dt
        if self._has_body:
            self._body.time = self._time
            self._bound.update_strengths_unsteady(dt, self._Uinfty, self._wake)
            self._wake.append(*self._bound.get_newly_shed())


class ExplicitEuler(Timestepper):
    """Timestepper using the explicit Euler method"""

    def advance(self, dt=None):
        """Advance the simulation for one timestep"""
        if not dt:
            dt = self._dt
        vel = self._wake_velocity()
        self._update_flow(self.wake.positions + vel * dt, dt)

class RungeKutta2(Timestepper):
    """Timestepper using 2nd-order Runge Kutta"""

    def advance(self, dt=None):
        """Advance the solution for one timestep"""
        if not dt:
            dt = self._dt
        x = self.wake.positions
        k1 = self._wake_velocity()
        k2 = self._wake_velocity(x + dt/2 * k1, dt/2)
        self._update_flow(x + dt * k2, dt)

class RungeKutta4(Timestepper):
    """Timestepper using 4th-order Runge Kutta"""

    def advance(self, dt=None):
        """Standard rk4"""
        if not dt:
            dt = self._dt
        x = self.wake.positions
        k1 = self._wake_velocity()
        k2 = self._wake_velocity(x + dt/2 * k1, dt/2)
        k3 = self._wake_velocity(x + dt/2 * k2, dt/2)
        k4 = self._wake_velocity(x + dt * k3, dt)
        self._update_flow(x + dt/6 * (k1 + 2 * k2 + 2 * k3 + k4), dt)
