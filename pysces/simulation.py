"""A module to easily set up and manage a simulation"""

__all__ = ['Simulation']

class Simulation(object):
    """Set up and run a boundary element simulation"""

    def __init__(self, body, Uinfty, dt, body_cls, wake_cls):
        """Initialize a simulation"""
        self._body = body
        self._Uinfty = Uinfty
        self._dt = dt
        self._body_panels = body_cls(body)
        self._wake = wake_cls()
        self.initialize()

    def initialize(self):
        """Initialize a simulation

        Solve for panel strengths so that surface boundary conditions are
        satisfied, and shed a particle into the wake so that overall circulation
        is zero.

        """
        self._time = 0
        self._body.time = 0
        # self._wake.reset()
        self._body_panels.update_strengths_unsteady(self._dt, self._Uinfty)
        self._wake.append(*self._body_panels.get_newly_shed())

    def advance(self, dt=None):
        """Advance the simulation for one timestep

        Use the explicit Euler method for advecting wake vortices

        """
        if not dt:
            dt = self._dt
        self._wake.advect(dt, self._Uinfty, self._body_panels.vortices)
        self._time += dt
        self._body.time = self._time
        self._body_panels.update_strengths_unsteady(dt, self._Uinfty, self._wake)
        self._wake.append(*self._body_panels.get_newly_shed())

    @property
    def time(self):
        """Current simulation time"""
        return self._time

    @property
    def body_panels(self):
        """Body panels used in the simulation"""
        return self._body_panels

    @property
    def wake_panels(self):
        """Wake panels used in the simulation"""
        return self._wake
