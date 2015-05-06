import unittest
from pysces.timestepper import *
from pysces.body import flat_plate
from pysces.panel import BoundVortices
from pysces.vortex import Vortices
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

class TestTimestepper(unittest.TestCase):
    def test_euler(self):
        body = flat_plate(20)
        Uinfty = (1,0)
        dt = 0.1
        flow = ExplicitEuler(dt, Uinfty, body, BoundVortices)
        self.assertEqual(flow.time, 0)
        self.assertEqual(len(flow.wake), 1)
        vort = flow.bound.vortices
        self.assertEqual(vort.circulation, -flow.wake.circulation)
        flow.advance()
        self.assertEqual(flow.time, dt)
        self.assertEqual(len(flow.wake), 2)

    def test_rk4(self):
        body = flat_plate(20)
        Uinfty = (1,0)
        dt = 0.1
        flow = RungeKutta4(dt, Uinfty, body, BoundVortices)
        self.assertEqual(flow.time, 0)
        self.assertEqual(len(flow.wake), 1)
        vort = flow.bound.vortices
        self.assertEqual(vort.circulation, -flow.wake.circulation)
        flow.advance()
        self.assertEqual(flow.time, dt)
        self.assertEqual(len(flow.wake), 2)

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
        self.assertLess(err, tol)

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
