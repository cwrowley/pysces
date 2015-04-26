import naca
import numpy as np

class Body(object):
    def __init__(self):
        self._time = 0

    @property
    def time(self):
        return self._time

    @time.setter
    def time(self, value):
        self._time = value

class Circle(Body):
    """Circle
    """
    def __init__(self, radius, num_points):
        """Return a circle with specified radius and number of points
        """
        self._radius = radius
        th = np.linspace(0, 2 * np.pi, num_points)
        self._x = radius * np.cos(th)
        self._y = radius * np.sin(th)

    def get_points(self):
        return self._x, self._y


class Airfoil(Body):
    """NACA 4-digit series airfoil
    """
    def __init__(self, code, num_points, **kwargs):
        """Return a NACA 4-digit series airfoil
        """
        self._x, self._y = naca.naca4(code, num_points, **kwargs)

    def get_points(self):
        return self._x, self._y


class BodyTransformation(Body):
    """Base class for transformations of existing bodies
    """
    def __init__(self, body):
        super(BodyTransformation, self).__init__()
        self._parent = body

    @property
    def time(self):
        Body.time.fget(self)

    @time.setter
    def time(self, value):
        Body.time.fset(self, value)
        # pass time to parent as well
        self._parent.time = value

    def get_points_from_parent(self):
        return self._parent.get_points()


class Translation(BodyTransformation):
    """Translation of an existing body
    """
    def __init__(self, body, displacement=(0,0)):
        super(Translation, self).__init__(body)
        self.displacement = displacement

    @property
    def displacement(self):
        return self._displacement

    @displacement.setter
    def displacement(self, value):
        self._displacement = np.array(value)

    def get_points(self):
        x, y = self.get_points_from_parent()
        q = np.vstack([x, y])
        if self.displacement.any():
            q += self.displacement[:, np.newaxis]
        return q[0,:], q[1,:]


class Rotation(BodyTransformation):
    """Rotation of an existing body
    """
    def __init__(self, body, angle=0.):
        super(Rotation, self).__init__(body)
        self.angle = angle

    @property
    def angle(self):
        return self._angle

    @angle.setter
    def angle(self, value):
        self._angle = value
        th = value * np.pi / 180
        self._R = np.array([[np.cos(th), np.sin(th)],[-np.sin(th), np.cos(th)]])

    def get_points(self):
        x, y = self.get_points_from_parent()
        q = np.vstack([x, y])
        if self._angle:
            q = np.dot(self._R, q)
        return q[0,:], q[1,:]


class Pitching(BodyTransformation):
    """Sinusoidal pitching for an existing body
    """
    def __init__(self, body, amplitude, frequency, phase=0.):
        super(Pitching, self).__init__(body)
        self._amplitude = amplitude
        self._frequency = frequency
        self._phase = phase * np.pi / 180
        self._body = Rotation(body)

    def get_points(self):
        angle = self._amplitude * np.sin(self._frequency * self._time
                                         + self._phase)
        self._body.angle = angle
        return self._body.get_points()


class Heaving(BodyTransformation):
    """Sinusoidal heaving for an existing body
    """
    def __init__(self, body, displacement, frequency, phase=0.):
        super(Heaving, self).__init__(body)
        self._displacement = np.array(displacement)
        self._frequency = frequency
        self._phase = phase * np.pi / 180
        self._body = Translation(body)

    def get_points(self):
        displacement = self._displacement * np.sin(self._frequency * self._time
                                                   + self._phase)
        self._body.displacement = displacement
        return self._body.get_points()


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
