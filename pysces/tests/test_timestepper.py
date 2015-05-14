import unittest
import sys
from pysces.timestepper import *
from pysces.body import flat_plate
from pysces.panel import BoundVortices
from pysces.vortex import Vortices
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

class TestTimestepper(unittest.TestCase):
    def check_timestepper(self, timestepper_cls):
        body = flat_plate(20)
        bound = BoundVortices(body)
        Uinfty = (1,0)
        dt = 0.1
        flow = timestepper_cls(dt, Uinfty, bound)
        self.assertEqual(flow.time, 0)
        self.assertEqual(len(flow.wake), 1)
        vort = flow.bound.vortices
        wake = flow.wake
        self.assertEqual(vort.circulation, -wake.circulation)
        flow.advance()
        self.assertEqual(flow.time, dt)
        self.assertEqual(len(wake), 2)
        self.assertEqual(vort.circulation, -wake.circulation)

    def test_euler(self):
        self.check_timestepper(ExplicitEuler)

    def test_rk2(self):
        self.check_timestepper(RungeKutta2)

    def test_rk4(self):
        self.check_timestepper(RungeKutta4)

    def check_vortex_pair(self, cls, tol):
        # compare with exact solution for a pair of vortices:
        # uniform rotation at frequency omega about center of vorticity (here 0)
        dt = 0.2
        Uinfty = (0,0)
        v1 = (-2,0)
        v2 = (1,0)
        g1 = 2 * np.pi
        g2 = 4 * np.pi
        vort = Vortices([v1,v2],[g1,g2])
        omega = 1./3
        stepper = cls(dt, Uinfty, wake=vort)
        num_steps = 10
        t = dt * num_steps
        v = np.array([np.cos(omega * t), np.sin(omega * t)])
        exact = np.array([-2 * v, v])
        for i in range(num_steps):
            stepper.advance()
        q = stepper.wake.positions
        err = np.sum((q - exact) * (q - exact))
        if sys.version_info >= (2, 7):
            self.assertLess(err, tol)
        else:
            # assertLess not available in Python 2.6
            self.assertTrue(err < tol)

    def test_vortex_pair_euler(self):
        tol = 0.07
        self.check_vortex_pair(ExplicitEuler, tol)

    def test_vortex_pair_rk2(self):
        tol = 1.6e-5
        self.check_vortex_pair(RungeKutta2, tol)

    def test_vortex_pair_rk4(self):
        tol = 6.4e-10
        self.check_vortex_pair(RungeKutta4, tol)

if __name__ == "__main__":
    unittest.main()
