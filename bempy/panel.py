from __future__ import division

import numpy as np

__all__ = ['BoundVortexPanels', 'FreeVortexParticles', 'SourceDoubletPanels']

class VortexPanels(object):
    pass

class BoundVortexPanels(object):
    """A class for bound vortex panels"""

    def __init__(self, body, Uinfty=(1,0)):
        self._body = body
        self._update(Uinfty)

    def _update(self, Uinfty=(1,0)):
        # here, Uinfty is used solely to determine direction of panel, for
        # placing colllocation points and vortex positions
        q = self._body.get_points(body_frame=True)
        dq = np.diff(q)
        self._numpanels = dq.shape[1]
        self._normals = (np.vstack([dq[1,:], -dq[0,:]]) /
                         np.linalg.norm(dq, axis=0))
        q25 = q[:,:-1] + 0.25 * dq
        q75 = q[:,:-1] + 0.75 * dq
        # vortex positions at 1/4 chord of panel
        # collocation points at 3/4 chord of panel
        # Determine orientation from Uinfty
        Uinfty = np.array(Uinfty)
        self._xvort = q25.copy()
        self._xcoll = q75.copy()
        # reverse direction where dq is against flow direction
        top, = np.where(np.dot(Uinfty, dq) <= 0)
        self._xvort[:,top] = q75[:,top]
        self._xcoll[:,top] = q25[:,top]
        self._gam = np.zeros(self._numpanels)
        self._influence_matrix = None

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
            for i, xv in enumerate(np.transpose(self._xvort)):
                vel = self.induced_velocity(self._xcoll, xv, 1)
                A[:, i] = np.sum(vel * self._normals, 0)
            self._influence_matrix = A
        return self._influence_matrix

    @staticmethod
    def induced_velocity(x, xvort, gam):
        r"""Compute velocity induced at points x by vortex at (xvort, gam)

        Induced velocity is

        .. math:: u_\theta = -\frac{\Gamma}{2 \pi r}
        """
        r = x - xvort[:,np.newaxis]
        rsq = np.sum(r * r, 0)
        return -gam / (2 * np.pi) * np.vstack([-r[1], r[0]]) / rsq

    def update_strengths(self, wake=None, Uinfty=(1,0)):
        rhs = self.compute_rhs(wake, Uinfty)
        self._gam = np.linalg.solve(self.influence_matrix, rhs)

    def compute_rhs(self, wake, Uinfty):
        # get collocation points and normals
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
            vel = np.zeros((2,self._numpanels))
        # assume body is not deforming: only motion is translation/rotation
        if motion:
            vel -= motion.map_velocity(self._xcoll)
        vel += np.array(Uinfty)[:,np.newaxis]
        # compute -v . n
        return -np.sum(vel * normals_inertial, 0)

    def add_wake_panel(self):
        pass

    def get_wake_panel(self):
        return None

    @property
    def vortices(self):
        return self._xvort, self._gam

    @property
    def collocation_pts(self):
        return self._xcoll

    @property
    def normals(self):
        return self._normals[0,:], self._normals[1,:]

class FreeVortexParticles(object):
    def __init__(self):
        pass

    def update(self, body, Uinfty, dt):
        pass

    def add_panels(self, panels):
        pass

    def induced_velocity(self, xcoll):
        return 0

    @property
    def vortices(self):
        return np.array([0,0]), 0


class SourceDoubletPanels(object):
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
