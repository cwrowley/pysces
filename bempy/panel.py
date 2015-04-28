from __future__ import division

import numpy as np

__all__ = ['BoundVortexPanels', 'FreeVortexParticles', 'SourceDoubletPanels']

class VortexPanels(object):
    pass

class BoundVortexPanels(object):
    def __init__(self, body):
        self._body = body
        q = body.get_points(body_frame=True)
        dq = np.diff(q)
        self._numpanels = dq.shape[1]
        self._normals = (np.vstack([dq[1,:], -dq[0,:]]) /
                         np.linalg.norm(dq, axis=0))
        q25 = q[:,:-1] + 0.25 * dq
        q75 = q[:,:-1] + 0.75 * dq
        # vortex positions at 1/4 chord of panel
        # collocation points at 3/4 chord of panel
        # assume first half goes from trailing edge to leading edge,
        #        second half from leading edge back to trailing edge
        half = self._numpanels // 2
        self._xvort = np.hstack([q75[:,:half], q25[:,half:]])
        self._xcoll = np.hstack([q25[:,:half], q75[:,half:]])
        self._gam = np.zeros(self._numpanels)

    def update_positions(self):
        self.panels = self._body.get_points()

    def update_strengths(self, wake, Uinfty, dt):
        # compute influence coefficients and RHS and solve for strengths
        pass

    def get_wake_panels(self):
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

    def get_wake_panels(self):
        return None
