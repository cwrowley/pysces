import numpy as np

class EuclideanTransformation(object):
    """A class representing elements of SE(2), the special Euclidean group"""

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
