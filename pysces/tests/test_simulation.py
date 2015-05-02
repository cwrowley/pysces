import unittest
from pysces.simulation import *
from pysces.body import flat_plate
from pysces.panel import BoundVortexPanels
from pysces.vortex import Vortices
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

class TestSimulation(unittest.TestCase):
    def test_simulation(self):
        body = flat_plate(20)
        Uinfty = (1,0)
        dt = 0.1
        flow = Simulation(body, Uinfty, dt, BoundVortexPanels, Vortices)
        self.assertEqual(flow.time, 0)
        self.assertEqual(len(flow.wake_panels), 1)
        vort = flow.body_panels.vortices
        self.assertEqual(vort.circulation, -flow.wake_panels.circulation)
        flow.advance()
        self.assertEqual(flow.time, dt)
        self.assertEqual(len(flow.wake_panels), 2)

if __name__ == "__main__":
    unittest.main()
