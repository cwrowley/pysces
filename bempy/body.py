import numpy as np
from euclid import EuclideanTransformation

__all__ = ['Body', 'Circle', 'Airfoil',
           'TransformedBody', 'Pitching', 'Heaving']

class Body(object):
    """Base class for representing bodies
    """
    def __init__(self):
        self._time = 0

    @property
    def time(self):
        return self._time

    @time.setter
    def time(self, value):
        self._time = value

    def get_body(self):
        """Return the Body object in the body-fixed frame
        """
        return self

    def get_transformation(self):
        """Return the transformation from the body-fixed to inertial frame
        """
        return None

class Circle(Body):
    """Circle
    """
    def __init__(self, radius, num_points):
        """Return a circle with specified radius and number of points
        """
        super(Circle, self).__init__()
        self._radius = radius
        th = np.linspace(0, 2 * np.pi, num_points)
        self._points = radius * np.vstack([np.cos(th), np.sin(th)])

    def get_points(self, **kwargs):
        return self._points


class Airfoil(Body):
    """NACA 4-digit series airfoil
    """
    def __init__(self, code, num_points, zero_thick_te=False, uniform=False):
        """Return a NACA 4-digit series airfoil
        """
        super(Airfoil, self).__init__()
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
        self._points = np.vstack([x, y])

    def get_points(self, **kwargs):
        return self._points


class TransformedBody(object):
    """Base class for Euclidean transformations of existing bodies
    """
    def __init__(self, body, angle=0, displacement=(0,0)):
        self._parent = body
        self._body = body.get_body()
        self._transformation = EuclideanTransformation(angle, displacement)

    def get_body(self):
        return self._body

    def get_transformation(self):
        return self._transformation.compose(self._parent.get_transformation())

    @property
    def displacement(self):
        return self._transformation.displacement

    @displacement.setter
    def displacement(self, value):
        self._transformation.displacement = value

    @property
    def angle(self):
        return self._transformation.angle

    @angle.setter
    def angle(self, value):
        self._transformation.angle = value

    @property
    def time(self):
        return self._body.time

    @time.setter
    def time(self, value):
        self._body.time = value

    def get_points(self, body_frame=False):
        q = self._body.get_points()
        if body_frame:
            return q
        return self.get_transformation().xform_position(q)


class Pitching(TransformedBody):
    """Sinusoidal pitching for an existing body
    """
    def __init__(self, body, amplitude, frequency, phase=0.):
        super(Pitching, self).__init__(body)
        self._amplitude = amplitude
        self._frequency = frequency
        self._phase = phase * np.pi / 180

    def get_transformation(self):
        self.angle = self._amplitude * np.sin(self._frequency * self.time
                                              + self._phase)
        return super(Pitching, self).get_transformation()


class Heaving(TransformedBody):
    """Sinusoidal heaving for an existing body
    """
    def __init__(self, body, displacement, frequency, phase=0.):
        super(Heaving, self).__init__(body)
        self._displacement = np.array(displacement)
        self._frequency = frequency
        self._phase = phase * np.pi / 180

    def get_transformation(self):
        displacement = self._displacement * np.sin(self._frequency * self.time
                                                   + self._phase)
        self.displacement = displacement
        return super(Heaving, self).get_transformation()
