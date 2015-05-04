"""A module to easily set up and manage a simulation"""
import numpy as np
from .vortex import Vortices

__all__ = ['ExplicitEuler']

class ExplicitEuler(object):
    """Set up and run a boundary element simulation"""

    def __init__(self, body, Uinfty, dt, body_cls):
        """Initialize a simulation"""
        self._body = body
        self._Uinfty = np.array(Uinfty)
        self._dt = dt
        self._bound = body_cls(body)
        self.initialize()

    def initialize(self):
        """Initialize a simulation

        Solve for panel strengths so that surface boundary conditions are
        satisfied, and shed a particle into the wake so that overall circulation
        is zero.

        """
        self._time = 0
        self._body.time = 0
        self._wake = Vortices()
        self._bound.update_strengths_unsteady(self._dt, self._Uinfty)
        self._wake.append(*self._bound.get_newly_shed())

    def advance(self, dt=None):
        """Advance the simulation for one timestep

        Use the explicit Euler method for advecting wake vortices

        """
        if not dt:
            dt = self._dt
        pos = self._wake.positions
        vel = self._wake.induced_velocity(pos)
        vel += self._bound.induced_velocity(pos)
        vel += self._Uinfty
        self._wake.positions += vel * dt
        self._time += dt
        self._body.time = self._time
        self._bound.update_strengths_unsteady(dt, self._Uinfty, self._wake)
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
