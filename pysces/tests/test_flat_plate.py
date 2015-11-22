import numpy as np
import unittest
from pysces import *

class TestFlatPlate(unittest.TestCase):
    """ Test unsteady flat plate problem for small angles of attack by 
    compareing with exact results from thin airfoil theory.  The L2 norm of the
    computed solution relative to the exact solution must be below a given
    threshold."""

    _threshold = 1.e-2

    def test_flat_plate(self):
        aoa_deg = [2,4,6]
        for aoa in aoa_deg:
            err_norm = self.flat_plate(aoa)
            success = True if (err_norm< self._threshold) else False
            self.assertTrue(success)

    def flat_plate(self, aoa_deg = 2.0):
        num_points = 20
        dt = .5
        Uinfty = (1,0)
        num_steps = 100

        aoa_rad = aoa_deg*np.pi/180
        body = TransformedBody(flat_plate(num_points), angle=aoa_deg)
        bound = BoundVortices(body)
        stepper = RungeKutta2(dt, Uinfty, bound)
        for i in range(num_steps):
            stepper.advance()
        vort = stepper.wake.positions
        s, dgam = self.compute_gam(body, stepper.bound.vortices)
        # Exact solution given on p.123 of Kuethe and Chow
        dgam_exact = 2 * aoa_rad * np.sqrt((1-s) / s)

        # Less accuracy at leading edge, so only compute along last half of
        # chord; since points advance from trailing edge to leading edge, we
        # want the first half of the array.
        n = len(s)//2
        err_mag = np.linalg.norm(dgam[0:n]-dgam_exact[0:n])
        return err_mag

    def compute_gam(self, body, vort):
        q = body.get_points()
        gam = vort.strengths
        ds = -np.linalg.norm(np.diff(q, axis=0), axis=1)
        dgam = gam / ds
        xvort = vort.positions
        s = np.sqrt(xvort[:,0]**2 + xvort[:,1]**2)
        return s, dgam

if __name__ == '__main__':
    unittest.main()
