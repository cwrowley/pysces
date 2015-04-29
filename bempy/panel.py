from __future__ import division

import numpy as np

__all__ = ['BoundVortexPanels', 'FreeVortexParticles', 'SourceDoubletPanels']

class VortexPanels(object):
    pass

class BoundVortexPanels(object):
    def __init__(self, body, Uinfty=(1,0)):
        # Uinfty is used solely to determine direction of panel, for placing
        # colllocation points and vortex positions
        self._body = body
        q = body.get_points(body_frame=False)
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

    def update_positions(self):
        self.panels = self._body.get_points()

    @staticmethod
    def induced_velocity(x, xvort, gam):
        """Compute velocity induced at points x by vortex at (xvort, gam)

        Induced velocity is u_theta = -gam / (2 pi r)"""
        r = x - xvort[:,np.newaxis]
        rsq = np.sum(r * r, 0)
        return -gam / (2 * np.pi) * np.vstack([-r[1], r[0]]) / rsq

    def update_strengths(self, wake=None, Uinfty=(1,0)):
        # compute influence coefficients and RHS and solve for strengths
        rhs = self.compute_rhs(wake, Uinfty)
        A = self.compute_influence_matrix()
        self._gam = np.linalg.solve(A, rhs)

    def compute_rhs(self, wake, Uinfty):
        # velocity induced by wake
        if wake:
            vel = wake.induced_velocity(self._xcoll)
        else:
            vel = np.zeros((2,self._numpanels))
        # assume body is not deforming: only motion is translation/rotation
        motion = self._body.get_transformation()
        if motion:
            vel -= motion.map_velocity(self._xcoll)
        vel += np.array(Uinfty)[:,np.newaxis]
        # compute -v . n
        return -np.sum(vel * self._normals, 0)

    def compute_influence_matrix(self):
        n = self._numpanels
        A = np.zeros((n, n))
        for i, xv in enumerate(np.transpose(self._xvort)):
            vel = self.induced_velocity(self._xcoll, xv, 1)
            A[:, i] = np.sum(vel * self._normals, 0)
        return A

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
