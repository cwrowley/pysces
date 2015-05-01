import numpy as np

__all__ = ['RigidMotion']

class RigidMotion(object):
    """A class representing rigid body motions, elements of TSE(2)"""

    def __init__(self, theta, x, thetadot=0, xdot=(0,0)):
        """angles are in radians"""
        self._theta = theta
        self._x = np.array(x)
        self._thetadot = thetadot
        self._xdot = np.array(xdot)
        self._update()

    def __repr__(self):
        if self._thetadot or self._xdot.any():
            return ("RigidMotion(%s, (%s, %s), %s, (%s, %s))" %
                    tuple(map(str, (self._theta, self._x[0], self._x[1],
                    self._thetadot, self._xdot[0], self._xdot[1]))))
        else:
            return ("RigidMotion(%s, (%s, %s))" %
                    tuple(map(str, (self._theta, self._x[0], self._x[1]))))

    def __str__(self):
        return self.__repr__()

    @property
    def theta(self):
        return self._theta

    @theta.setter
    def theta(self, value):
        self._theta = value
        self._update()

    def _update(self):
        c = np.cos(self._theta)
        s = np.sin(self._theta)
        self._R = np.array([[c, -s], [s, c]])
        self._Rdot = np.array([[-s, -c], [c, -s]]) * self._thetadot

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, value):
        self._x = np.array(value)

    @property
    def thetadot(self):
        return self._thetadot

    @property
    def xdot(self):
        return self._xdot

    @classmethod
    def identity(cls):
        return cls(0, (0,0))

    def __eq__(self, other):
        return (np.array_equal(self._x, other._x) and
                np.array_equal(self._xdot, other._xdot) and
                self._theta == other._theta and
                self._thetadot == other._thetadot)

    def __ne__(self, other):
        return not self == other

    def inverse(self):
        """Return the inverse element"""
        xinv = -np.dot(self._R.T, self._x)
        xdotinv = -np.dot(self._Rdot.T, self._x) - np.dot(self._R.T, self._xdot)
        return RigidMotion(-self._theta, xinv, -self._thetadot, xdotinv)

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
        theta = self._theta + other._theta
        x = np.dot(self._R, other._x) + self._x
        thetadot = self._thetadot + other._thetadot
        xdot = (np.dot(self._Rdot, other._x) + np.dot(self._R, other._xdot) +
                self._xdot)
        return RigidMotion(theta, x, thetadot, xdot)

    def map_position(self, q):
        """Return the action of the transformation on the given vector(s) q

        `q` can be an array of length 2 or a 2d-array with 2 rows

        The group action of the element (R, x) is given by
            q -> R . q + x
        """
        if self._theta:
            q_new = np.dot(q, np.transpose(self._R))
        else:
            q_new = np.array(q, copy=True)
        if self._x.any():
            q_new += self._x
            # if q.ndim == 1:
                # q_new += self._x
            # else:
                # q_new += self._x[:, np.newaxis]
        return q_new

    def map_vector(self, qdot):
        """Return the action of the transformation on the tangent vector(s) qdot

        `qdot` can be an array of length 2 or a 2d-array with 2 rows

        The tangent action is given by
            qdot -> R . qdot
        """
        if self._theta:
            return np.dot(qdot, np.transpose(self._R))
        else:
            return np.array(qdot, copy=True)

    def map_velocity(self, q, qdot=None):
        """Return the velocity of the transformed base point q, velocity qdot

        If transformation is (R, x), then
            d/dt (R, x) q = Rdot q + R qdot + xdot
        """
        qdot_new = np.zeros_like(q)
        if self._thetadot:
            qdot_new += np.dot(q, np.transpose(self._Rdot))
        if self._theta and qdot is not None and qdot.any():
            qdot_new += np.dot(qdot, np.transpose(self._R))
        if self._xdot.any():
            qdot_new += self._xdot
            # if q.ndim == 1:
                # qdot_new += self._xdot
            # else:
                # qdot_new += self._xdot[:, np.newaxis]
        return qdot_new
