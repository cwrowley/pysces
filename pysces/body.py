import numpy as np
from .motion import RigidMotion

__all__ = ['Body', 'TransformedBody', 'Pitching', 'Heaving',
           'cylinder', 'flat_plate', 'naca_airfoil']

class Body(object):
    """Base class for representing bodies
    """
    def __init__(self, points):
        """Create a body with nodes at the given points

        Parameters
        ----------
        points : 2d array, shape (n,2)
            Array of points defining the boundary of the body
            For a closed body, the boundary curve should be positively oriented
            (counter-clockwise around outside of body), starting from trailing
            edge
        """
        self._time = 0
        self._points = points

    @property
    def time(self):
        """The time used to specify the body's motion"""
        return self._time

    @time.setter
    def time(self, value):
        self._time = value

    def get_points(self, body_frame=False):
        return self._points

    def get_body(self):
        """Return the Body object in the body-fixed frame"""
        return self

    def get_motion(self):
        """Return the transformation from the body-fixed to inertial frame"""
        return None

def cylinder(radius, num_points):
    """Return a circular Body with the given radius and number of points"""
    th = np.linspace(0, 2 * np.pi, num_points)
    points = radius * np.array([np.cos(th), np.sin(th)]).T
    return Body(points)

def flat_plate(num_points):
    """Return a flat plate with the given number of points"""
    x = np.linspace(1, 0, num_points)
    y = np.zeros_like(x)
    return Body(np.array([x, y]).T)

def naca_airfoil(code, num_points, zero_thick_te=False, uniform=False):
    """Return a NACA 4-digit series airfoil"""
    # extract parameters from 4-digit code
    code_str = "%04d" % int(code)
    if len(code_str) != 4:
        raise ValueError("NACA designation is more than 4 digits")
    max_camber = 0.01 * int(code_str[0])
    p = 0.1 * int(code_str[1])  # location of max camber
    thickness = 0.01 * int(code_str[2:])
    if uniform:
        x = np.linspace(0, 1, num_points)
    else:
        # closer spacing near leading edge
        theta = np.linspace(0, 0.5 * np.pi, num_points)
        x = 1 - np.cos(theta)

    # thickness
    coefs = [-0.1015, 0.2843, -0.3516, -0.1260, 0, 0.2969]
    if zero_thick_te:
        coefs[0] = -0.1036
    y_thick = 5 * thickness * (np.polyval(coefs[:5], x) +
                               coefs[5] * np.sqrt(x))

    # camber
    front = np.where(x <= p)
    back = np.where(x > p)
    y_camber = np.zeros_like(x)
    if p:
        y_camber[front] = max_camber * x[front] / p**2 * (2 * p - x[front])
        y_camber[back] = max_camber * ((1. - x[back])/(1. - p)**2 *
                                       (1 + x[back] - 2 * p))
    x = np.hstack([x[-1:0:-1], x])
    y = np.hstack([y_camber[-1:0:-1] + y_thick[-1:0:-1],
                   y_camber - y_thick])
    return Body(np.array([x, y]).T)


class TransformedBody(object):
    """Base class for rigid (Euclidean) transformations of existing bodies
    """
    def __init__(self, body, angle=0, displacement=(0,0)):
        """angles are clockwise, in degrees"""
        self._parent = body
        self._body = body.get_body()
        self._motion = RigidMotion(-angle * np.pi / 180, displacement)

    def get_body(self):
        return self._body

    def get_motion(self):
        self._update()
        return self._motion.compose(self._parent.get_motion())

    def set_motion(self, value):
        self._motion = value

    @property
    def time(self):
        return self._body.time

    @time.setter
    def time(self, value):
        self._body.time = value

    def _update(self):
        # update body motion: subclasses override this
        pass

    def get_points(self, body_frame=False):
        q = self._body.get_points()
        if body_frame:
            return q
        return self.get_motion().map_position(q)


class Pitching(TransformedBody):
    """Sinusoidal pitching for an existing body
    """
    def __init__(self, body, amplitude, frequency, phase=0.):
        """amplitude and phase given in degrees"""
        super(Pitching, self).__init__(body)
        self._amplitude = amplitude * np.pi / 180
        self._frequency = frequency
        self._phase = phase * np.pi / 180

    def _update(self):
        theta = self._frequency * self.time + self._phase
        alpha = self._amplitude * np.sin(theta)
        alphadot = self._amplitude * self._frequency * np.cos(theta)
        self.set_motion(RigidMotion(-alpha, (0,0), -alphadot, (0,0)))


class Heaving(TransformedBody):
    """Sinusoidal heaving for an existing body
    """
    def __init__(self, body, displacement, frequency, phase=0.):
        super(Heaving, self).__init__(body)
        self._displacement = np.array(displacement)
        self._frequency = frequency
        self._phase = phase * np.pi / 180

    def _update(self):
        theta = self._frequency * self.time + self._phase
        x = self._displacement * np.sin(theta)
        xdot = self._displacement * self._frequency * np.cos(theta)
        self.set_motion(RigidMotion(0, x, 0, xdot))
