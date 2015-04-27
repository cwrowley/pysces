import naca
import numpy as np

class EuclideanTransformation(object):
    def __init__(self, angle, displacement):
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

    def inverse(self):
        return EuclideanTransformation(-self._angle,
                                       -np.dot(self._R.T, self._displacement))

    def compose(self, other):
        """Return the composition of self (left) with other (right)

        If g1 = (R1, v1)
           g2 = (R2, v2)
        where R1, R2 are rotation matrices and v1, v2 are displacement vectors,
        then
            g1 g2 = (R1 R2, R1 v2 + v1)
        """
        if other is None:
            return self
        angle = self.angle + other.angle
        displacement = np.dot(self._R, other.displacement) + self._displacement
        return EuclideanTransformation(angle, displacement)

    def xform_position(self, vector):
        """Return the action of the transformation on the given vectors

        `vector` can be an array of length 2 or a 2d-array with 2 rows
        """
        newvector = np.array(vector, copy=True)
        if self._angle:
            newvector = np.dot(self._R, newvector)
        if self._displacement.any():
            if vector.ndim == 1:
                newvector += self._displacement
            else:
                newvector += self._displacement[:, np.newaxis]
        return newvector


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
        super(Airfoil, self).__init__()
        self._x, self._y = naca.naca4(code, num_points, **kwargs)

    def get_points(self):
        return self._x, self._y


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

    def get_points_body_frame(self):
        return self._body.get_points()

    def get_points(self):
        x, y = self.get_points_body_frame()
        q = self.get_transformation().xform_position(np.vstack([x, y]))
        return q[0,:], q[1,:]


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


class VortexPanels(object):
    pass

class BoundVortexPanels(object):
    def __init__(self, body):
        self._body = body
        self.panels = body.get_points()

    def update_positions(self):
        self.panels = self._body.get_points()

    def update_strengths(self, wake, Uinfty, dt):
        # compute influence coefficients and RHS and solve for strengths
        pass

    def get_wake_panels(self):
        return None

class FreeVortexParticles(object):
    def __init__(self):
        pass

    def update(self, body, Uinfty, dt):
        pass

    def add_panels(self, panels):
        pass

    @property
    def vortices(self):
        return 0, 0, 0


class SourceDoubletPanels(object):
    def __init__(self, body):
        self._body = body
        self.panels = body.get_points()

    def update_positions(self):
        self.panels = self._body.get_points()

    def update_strengths(self, wake, Uinfty, dt):
        # compute influence coefficients and RHS and solve for strengths
        pass

    def get_wake_panels(self):
        return None


def time_advance(body, wake, Uinfty, dt):
    # todo: how does dt come in to updating strengths of body panels, and in
    # particular wake panel to be shed?

    # 1) define the wake panel: need Uinf and dt to determine length
    # 2) to determine rhs for update_strengths, need Uinf for relative vel
    # 3) wake update needs Uinf and dt (relative vel)
    body.update_strengths(wake, Uinfty, dt) # might need dt to determine length of TE panel
    wake.update(body, Uinfty, dt)
    shed_panels = body.get_wake_panels()
    wake.add_panels(shed_panels)

def compute_forces(body, wake):
    return 0, 0
