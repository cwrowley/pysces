import naca
import numpy as np

class Body(object):
    def __init__(self):
        pass

class Airfoil(Body):
    """NACA 4-digit series airfoil
    """
    def __init__(self, code, num_points, **kwargs):
        """Return a NACA 4-digit series airfoil
        """
        self._x, self._y = naca.naca4(code, num_points, **kwargs)

    def get_points(self):
        return self._x, self._y


class Translation(Body):
    """Translation/rotation of an existing body
    """
    def __init__(self, body, angle=0., displacement=None):
        self.body = body
        self.angle = angle
        self.displacement = displacement

    @property
    def angle(self):
        return self._angle

    @angle.setter
    def angle(self, value):
        self._angle = value
        th = value * np.pi / 180
        self._R = np.array([[np.cos(th), np.sin(th)],[-np.sin(th), np.cos(th)]])

    @property
    def displacement(self):
        return self._displacement

    @displacement.setter
    def displacement(self, value):
        self._displacement = np.array(value)

    def get_points(self):
        x, y = self.body.get_points()
        q = np.vstack([x, y])
        if self._angle:
            q = np.dot(self._R, q)
        if self.displacement.any():
            q += self.displacement[:, np.newaxis]
        return q[0,:], q[1,:]


class VortexPanels(object):
    pass

class BoundVortexPanels(VortexPanels):
    def __init__(self, body):
        self.panels = body.get_points()

    def update_positions(self, body):
        self.panels = body.get_points()

    def update_strengths(self, wake):
        # compute influence coefficients and RHS and solve for strengths
        pass

    def get_wake_panels(self):
        return None

class FreeVortexPanels(VortexPanels):
    def __init__(self):
        pass

    def update(self, dt):
        pass

    def add_panels(self, panels):
        pass

    @property
    def vortices(self):
        return 0, 0, 0


def time_advance(body, wake, dt):
    # todo: how does dt come in to updating strengths of body panels, and in
    # particular wake panel to be shed?
    body.update_strengths(wake)
    wake.update(dt)
    shed_panels = body.get_wake_panels()
    wake.add_panels(shed_panels)

def compute_forces(body, wake):
    return 0, 0
