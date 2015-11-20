from __future__ import division

import numpy as np
import sys
import matplotlib.pyplot as plt
from .vortex import Vortices

__all__ = ['BoundVortices', 'BoundSourceDoublets']

class BoundVortices(object):
    """A class for bound vortex panels"""

    def __init__(self, body, Uinfty=(1,0)):
        self._body = body
        self._time = 0
        self._update(Uinfty)

    def _update(self, Uinfty=(1,0)):
        # Uinfty is used here solely to determine direction of panel and for
        # distributing bound vortices and collocation points
        q = self._body.get_points(body_frame=True)
        dq = np.diff(q, axis=0)
        self._numpanels = dq.shape[0]
        self._tangents = dq / np.linalg.norm(dq, axis=1)[:,np.newaxis]
        self._normals = np.transpose(np.array([dq[:,1], -dq[:,0]]) /
                                     np.linalg.norm(dq, axis=1))

        q25 = q[:-1] + 0.25 * dq
        q75 = q[:-1] + 0.75 * dq
        # vortex positions at 1/4 chord of panel
        # collocation points at 3/4 chord of panel
        # Determine orientation from Uinfty
        Uinfty = np.array(Uinfty)
        xvort = q25.copy()
        self._xcoll = q75.copy()
        # reverse direction where dq is against flow direction
        top, = np.where(np.dot(dq, Uinfty) <= 0)
        xvort[top] = q75[top]
        self._xcoll[top] = q25[top]
        # find trailing edge and wake vortex direction
        if np.linalg.norm(q[0] - q[-1]) < 0.005:
            # closed body
            self._trailing_edge = 0.5 * (q[0] + q[-1])
            wake_dir = dq[-1] - dq[0]
            self._wake_dir = wake_dir / np.linalg.norm(wake_dir)
        else:
            # thin airfoil
            self._trailing_edge = q[0]
            self._wake_dir = -dq[0] / np.linalg.norm(dq[0])
        self._vortices = Vortices(xvort)
        self._influence_matrix = None

        # Uncomment the following lines to see the vortex and collocation
        # points.
        #plt.plot(q[:,0],q[:,1],'-k')
        #plt.plot(xvort[:,0],xvort[:,1],'ro')
        #plt.plot(self._xcoll[:,0],self._xcoll[:,1],'bs')
        #plt.axis('equal')
        #plt.legend(['body','vortices','collocation pts'])
        #plt.show()

    def update_positions(self):
        # If non-rigid bodies are used, update panel positions here.
        #
        # Note that if the only motion is rigid body motion, the panel positions
        # do not need to be updated, since they are in body-fixed frame

        # need to recompute influence matrix when points change
        self._influence_matrix = None

    @property
    def influence_matrix(self):
        if self._influence_matrix is None:
            # time to recompute
            n = self._numpanels
            A = np.zeros((n, n))
            for j, vort in enumerate(self._vortices):
                vel = self._vortices.induced_velocity_single(self._xcoll,
                                                             vort[0], 1)
                A[:, j] = np.sum(vel * self._normals, 1)
            self._influence_matrix = A
        return self._influence_matrix

    @property
    def num_panels(self):
        return self._numpanels

    @property
    def tangents(self):
        return self._tangents

    @property
    def normals(self):
        return self._normals

    def update_strengths(self, Uinfty=(1,0)):
        """Update vortex strengths"""
        rhs = self.compute_rhs(Uinfty)
        self._vortices.strengths = np.linalg.solve(self.influence_matrix, rhs)

    def update_strengths_unsteady(self, dt, Uinfty=(1,0), wake=None, circ=None,
                                  wake_fac=0.25):
        """Update strengths for unsteady calculation

        Shed a new wake panel (not added into wake)

        Parameters
        ----------
        dt : float
            Timestep
        Uinfty : array_like, optional
            Farfield fluid velocity (default (1,0))
        wake : wake panel object, optional
            Induces velocities on the body (default None)
        circ : float, optional
            Total bound circulation, for enforcing Kelvin's circulation theorem.
            If None (default), obtain the total circulation from the wake,
            assuming overall circulation (body + wake) is zero
        wake_fac : float, optional
            New wake vortex is placed a distance wake_fac * Uinfty * dt from
            trailing edge (see Katz & Plotkin, p390).
        """

        # determine new wake vortex position (in body-fixed frame)
        distance = wake_fac * np.sqrt(Uinfty[0]**2 + Uinfty[1]**2) * dt
        x_shed = self._trailing_edge + distance * self._wake_dir
        # compute velocity induced on collocation points by newly shed vortex
        # (done in the body-fixed frame)
        shed_vel = self._vortices.induced_velocity_single(self._xcoll, x_shed, 1)
        shed_normal = np.sum(shed_vel * self._normals, 1)
        # determine overall influence matrix, including newly shed vortex
        # last equation: sum of all the vortex strengths = total circulation
        A = np.vstack([np.hstack([self.influence_matrix,
                                  shed_normal[:,np.newaxis]]),
                       np.ones((1, self._numpanels + 1))])

        rhs0 = self.compute_rhs(Uinfty, wake)
        if circ is None:
            if wake is None:
                circ = 0
            else:
                circ = -wake.circulation
        rhs = np.hstack([rhs0, circ])

        gam = np.linalg.solve(A, rhs)
        self._vortices.strengths = gam[:-1]
        self._x_shed = x_shed
        self._gam_shed = gam[-1]

    def compute_rhs(self, Uinfty=(1,0), wake=None):
        # get collocation points and normals
        # if a motion is present, use it to map the collocation points and 
        # their normals from the body frame to the inertial frame.
        motion = self._body.get_motion()
        if motion:
            xcoll_inertial = motion.map_position(self._xcoll)
            normals_inertial = motion.map_vector(self._normals)
        else:
            xcoll_inertial = self._xcoll
            normals_inertial = self._normals
        # velocity induced by wake
        if wake:
            vel = wake.induced_velocity(xcoll_inertial)
        else:
            vel = np.zeros((self._numpanels, 2))
        # assume body is not deforming: only motion is translation/rotation
        if motion:
            vel -= motion.map_velocity(self._xcoll)
        vel += np.array(Uinfty)
        # compute -v . n
        return -np.sum(vel * normals_inertial, 1)

    def get_newly_shed(self):
        """Return newly shed wake vortex in the inertial frame

        Returns
        -------
        x_shed : 1d array, shape (2,)
            Location of newly shed wake vortex, in inertial frame
        gam_shed : float
            Strength of newly shed vortex
        """
        motion = self._body.get_motion()
        if motion:
            x_shed_inertial = motion.map_position(self._x_shed)
        else:
            x_shed_inertial = np.array(self._x_shed, copy=True)
        return x_shed_inertial, self._gam_shed

    def induced_velocity(self, x):
        return self._vortices.induced_velocity(x, self._body.get_motion())

    @property
    def time(self):
        return self._time

    @time.setter
    def time(self, value):
        self._time = value
        self._body.time = value

    @property
    def vortices(self):
        return self._vortices

    @property
    def collocation_pts(self):
        return self._xcoll

    @property
    def normals(self):
        # return self._normals[0,:], self._normals[1,:]
        return self._normals


class BoundSourceDoublets(object):
    def __init__(self, body):
        self._body = body
        self.panels = body.get_points()

    def update_positions(self):
        self.panels = self._body.get_points()

    def update_strengths(self, wake, Uinfty, dt):
        # compute influence coefficients and RHS and solve for strengths
        pass

    def get_wake_panel(self):
        return None
