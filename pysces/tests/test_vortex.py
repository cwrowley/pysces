import unittest
from pysces.vortex import *
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
from pysces.motion import RigidMotion

class TestVortex(unittest.TestCase):
    def test_init_empty(self):
        vort = Vortices()
        self.assertEqual(len(vort), 0)

    def check_vortices(self, vort, pos, gam=None):
        n = len(pos)
        self.assertEqual(len(vort), n)
        for i in range(n):
            assert_array_equal(vort.positions[i], pos[i])
            if gam is None:
                self.assertEqual(vort.strengths[i], 0)
            else:
                self.assertEqual(vort.strengths[i], gam[i])

    def test_init_tuple(self):
        v1 = (0,0)
        vort = Vortices(v1)
        self.check_vortices(vort, [v1])

    def test_init_array1(self):
        v1 = np.array((0,0))
        vort = Vortices(v1)
        self.check_vortices(vort, [v1])

    def test_init_array2(self):
        v1 = (-1,0)
        v2 = (1,0)
        vort = Vortices(np.array([v1, v2]))
        self.check_vortices(vort, [v1, v2])

    def test_init_list(self):
        v1 = (-1,0)
        v2 = (1,0)
        vort = Vortices([v1, v2])
        self.check_vortices(vort, [v1, v2])

    def test_init_strength(self):
        v1 = (0,0)
        s1 = 1
        vort = Vortices([v1], [s1])
        self.check_vortices(vort, [v1], [s1])

    def test_init_strength_list(self):
        v1 = (-1,0)
        v2 = (1,0)
        g1 = 2
        g2 = 3
        pos = [v1, v2]
        gam = [g1, g2]
        vort = Vortices(pos, gam)
        self.check_vortices(vort, pos, gam)

    def test_append(self):
        v1 = (-1,0)
        v2 = (1,0)
        g1 = 2
        g2 = 3
        vort = Vortices([v1], [g1])
        self.check_vortices(vort, [v1], [g1])
        vort.append(v2, g2)
        self.check_vortices(vort, [v1,v2], [g1,g2])

    def test_append_empty(self):
        vort = Vortices()
        v = (0,0)
        s = 1
        vort.append(v, s)
        self.check_vortices(vort, [v], [s])

    def test_iter(self):
        v1 = (-1,0)
        v2 = (1,0)
        pos = [v1, v2]
        gam = [1,2]
        vort = Vortices(pos, gam)
        for i, v in enumerate(vort):
            assert_array_equal(v[0], pos[i])
            self.assertEqual(v[1], gam[i])

    def test_iter_empty(self):
        for x, gam in Vortices():
            pass

    def test_circulation(self):
        pos = [(0,0), (1,0)]
        gam = [42, -13]
        vort = Vortices(pos, gam)
        self.assertEqual(vort.circulation, 42 - 13)

    def test_core_radius(self):
        vort = Vortices()
        default_rad = 1.e-3
        self.assertEqual(vort.core_radius, default_rad)
        vort.core_radius = 0.5
        self.assertEqual(vort.core_radius, 0.5)
        Vortices.core_radius = 0.2
        vort2 = Vortices()
        self.assertEqual(vort2.core_radius, 0.2)

    def test_induced_velocity_single_tuple(self):
        v = (2,0)
        x = (3,0)
        gam = 2 * np.pi
        vel = Vortices().induced_velocity_single(x, v, gam)
        assert_array_equal(vel, (0,-1))

    def test_induced_velocity_single_array(self):
        v = (0,0)
        x = np.array([(1,0), (2,0), (3,0)])
        gam = 2 * np.pi
        vel = Vortices().induced_velocity_single(x, v, gam)
        vel_expected = np.array([(0, -1), (0, -1./2), (0, -1./3)])
        assert_array_equal(vel, vel_expected)

    def test_induced_velocity_tuple(self):
        v1 = (-1,0)
        v2 = (1,0)
        g1 = -2 * np.pi
        g2 = 2 * np.pi
        vort = Vortices([v1, v2], [g1, g2])
        x = (0,0)
        vel = vort.induced_velocity(x)
        vel_expected = (0,2)
        assert_array_equal(vel, vel_expected)

    def test_induced_velocity_array(self):
        v1 = (-1,0)
        v2 = (1,0)
        g1 = -2 * np.pi
        g2 = 2 * np.pi
        vort = Vortices([v1, v2], [g1, g2])
        x = np.array([(0,-1), (0,0), (0,1)])
        vel = vort.induced_velocity(x)
        vel_expected = np.array([(0,1), (0,2), (0,1)])
        assert_array_equal(vel, vel_expected)

    def test_induced_velocity_motion(self):
        v = (1,0)
        vort = Vortices(v, 2 * np.pi)
        motion = RigidMotion(np.pi/2, (0,0))
        x = (0,0)
        vel = vort.induced_velocity(x, motion=motion)
        vel_expected = (-1,0)
        assert_array_almost_equal(vel, vel_expected)

    def test_regularization(self):
        eps = 1.e-2
        vort = Vortices((0,0), 2 * np.pi)
        vort.core_radius = eps
        x0 = np.array((eps,0))
        x = np.array([0.5 * x0, x0, 2 * x0])
        vel = vort.induced_velocity(x)
        vel_expected = np.array([(0,-0.5/eps), (0,-1./eps), (0,-0.5/eps)])
        assert_array_equal(vel, vel_expected)
