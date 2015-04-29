from bempy.motion import *
import unittest
import numpy as np
# from numpy.testing import assert_array_equal, assert_array_almost_equal

class TestRigidMotion(unittest.TestCase):
    def assert_motion_almost_equal(self, g, h):
        np.testing.assert_almost_equal(g.theta, h.theta)
        np.testing.assert_array_almost_equal(g.x, h.x)
        np.testing.assert_almost_equal(g.thetadot, h.thetadot)
        np.testing.assert_array_almost_equal(g.xdot, h.xdot)

    def test_identity_action(self):
        e = RigidMotion(0, (0,0))
        q = (13,42)
        np.testing.assert_array_equal(q, e.map_position(q))

    def test_equality(self):
        g1 = RigidMotion(1, (2,3))
        g2 = RigidMotion(1, (2,3))
        self.assertEqual(g1, g2)
        g3 = RigidMotion(1, (2,3), thetadot=1)
        self.assertNotEqual(g1, g3)
        g4 = RigidMotion(1, (2,3), xdot=(1,0))
        self.assertNotEqual(g1, g4)
        g5 = RigidMotion(0, (2,3))
        self.assertNotEqual(g1, g5)
        g6 = RigidMotion(1, (2,4))
        self.assertNotEqual(g1, g6)

    def test_identity(self):
        e = RigidMotion(0, (0,0))
        self.assertEqual(e, RigidMotion.identity())

    def test_compose(self):
        g1 = RigidMotion(np.pi/2, (0,0))
        g2 = RigidMotion(0, (2,3))
        g3 = RigidMotion(np.pi/2, (2,3))
        self.assertEqual(g2.compose(g1), g3)
        g4 = RigidMotion(np.pi/2, (3,-2))
        self.assert_motion_almost_equal(g1.compose(g2), g4)

    def test_inverse(self):
        e = RigidMotion.identity()
        g1 = RigidMotion(1, (0,0))
        g1inv = RigidMotion(-1, (0,0))
        self.assertEqual(g1.inverse(), g1inv)

        g2 = RigidMotion(0, (2,3))
        g2inv = RigidMotion(0, (-2,-3))
        self.assertEqual(g2.inverse(), g2inv)

        g3 = RigidMotion(1,(2,3))
        g3inv = g1inv.compose(g2inv)
        self.assertEqual(g3.inverse(), g3inv)

    def test_map_position(self):
        rot = RigidMotion(np.pi/2, (1,2))
        x = np.array((13, 42))
        y = (43, -11)
        np.testing.assert_array_almost_equal(rot.map_position(x), y)

        xarr = np.array([[13,13],[42,42]])
        yarr = np.array([[43,43],[-11,-11]])
        np.testing.assert_array_almost_equal(rot.map_position(xarr), yarr)


    def test_map_velocity_inplace(self):
        rot = RigidMotion(0, (0,0), 3, (0,0))
        x = np.array((1, 0))
        v = np.array((0, -3))
        np.testing.assert_array_almost_equal(rot.map_velocity(x), v)
        np.testing.assert_array_almost_equal(rot.map_velocity(2*x), 2*v)

        vel = np.array((4,5))
        rot2 = RigidMotion(0, (0,0), 3, vel)
        np.testing.assert_array_almost_equal(rot2.map_velocity(x), v + vel)
        np.testing.assert_array_almost_equal(rot2.map_velocity(2*x), 2*v + vel)

    def test_map_velocity_rot(self):
        rot = RigidMotion(np.pi/2, (42,13), 3, (0,0))
        x = np.array((1, 0))
        v = np.array((-3, 0))
        np.testing.assert_array_almost_equal(rot.map_velocity(x), v)
        np.testing.assert_array_almost_equal(rot.map_velocity(2*x), 2*v)

        vel = np.array((4,5))
        rot2 = RigidMotion(np.pi/2, (42,13), 3, vel)
        np.testing.assert_array_almost_equal(rot2.map_velocity(x), v + vel)
        np.testing.assert_array_almost_equal(rot2.map_velocity(2*x), 2*v + vel)
