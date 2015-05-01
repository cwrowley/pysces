import numpy as np

__all__ = ['Vortices']

class Vortices(object):
    core_radius = 1.e-3

    def __init__(self, positions=None, strengths=None):
        if positions is None:
            self._positions = None
        else:
            self._positions = np.array(positions, ndmin=2, dtype=np.float64)

        self._circulation = 0
        if strengths is None:
            if positions is None:
                self._strengths = None
            else:
                self._strengths = np.zeros(self._positions.shape[0],
                                           dtype=np.float64)
        else:
            self._strengths = np.array(strengths, ndmin=1, dtype=np.float64)
            self._circulation = np.sum(self._strengths)

    @property
    def positions(self):
        return self._positions

    @property
    def strengths(self):
        return self._strengths

    @strengths.setter
    def strengths(self, value):
        self._strengths = np.array(value, ndmin=1, dtype=np.float64)

    @property
    def circulation(self):
        return self._circulation

    def __len__(self):
        if self._positions is None:
            return 0
        return self._positions.shape[0]

    def __iter__(self):
        if self._positions is None:
            return iter([])
        return iter(zip(self._positions, self._strengths))

    def append(self, position, strength):
        position = np.array(position, ndmin=2)
        strength = np.array(strength, ndmin=1)
        if self._positions is None:
            self._positions = position
            self._strengths = strength
            self._circulation = np.sum(strength)
        else:
            self._positions = np.append(self._positions, position, axis=0)
            self._strengths = np.append(self._strengths, strength)
            self._circulation += np.sum(strength)

    def induced_velocity_single(self, x, xvort, gam):
        r"""Compute velocity induced at points x by a single vortex

        Parameters
        ----------
        x : 2d array
            Locations at which to compute induced velocity.  Expressed as
            column vectors (i.e., shape should be (n,2))
        xvort : 1d array
            Location of vortex (shape should be (2,))
        gam : float
            Strength of vortex

        Notes
        -----
        Induced velocity is

        .. math:: u_\theta = -\frac{\Gamma}{2 \pi r}

        where r is the distance between the point and the vortex.  If this
        distance is less than :class:`core_radius` :math:`r_0`, the velocity is
        regularized as solid-body rotation, with

        .. math:: u_\theta = -\frac{\Gamma r}{2\pi r_0^2}
        """
        r = np.array(x, ndmin=2) - np.array(xvort)
        rsq = np.maximum(np.sum(r * r, 1), self.core_radius**2)
        # alternative regularization (Krasny, Eldredge)
        # rsq = np.sum(r * r, 1) + self.core_radius**2
        vel = np.array(r[:,[1,0]], copy=True)
        vel = gam / (2 * np.pi) * vel / rsq[:,np.newaxis]
        vel[:,1] *= -1
        # print(vel)
        return np.squeeze(vel)

    def induced_velocity(self, x, motion=None):
        """Compute the induced velocity at the given point(s)"""
        if motion is None:
            positions = self._positions
        else:
            positions = motion.map_position(self._positions)
        x = np.array(x)
        vel = np.zeros_like(x, dtype=np.float64)
        for xvort, gam in zip(positions, self._strengths):
            vel += self.induced_velocity_single(x, xvort, gam)
        return vel

    def advect(self, dt, Uinfty=(0,0), other=None):
        """Advect the vortex particles forward one step in time

        Parameters
        ----------
        dt : float
            Timestep
        Uinfty : array_like, optional
            Farfield velocity, default (0,0)
        ither : Vortex object (optional)
            Optional body also contributing to induced velocity of the particles

        Notes
        -----
        An explicit Euler update is used, where the particle positions are
        incremented by ``vel * dt``, where ``vel`` is the induced velocity

        """
        # explicit Euler update
        vel = self.induced_velocity(self._positions)
        if other:
            vel += other.induced_velocity(self._positions)
        Uinfty = np.array(Uinfty)
        if Uinfty.any():
            vel += Uinfty
        self._positions += vel * dt

    # def reset(self):
    #     self._positions = None
    #     self._strengths = None
    #     self.
