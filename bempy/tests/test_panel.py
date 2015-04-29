import unittest
from bempy.body import Body
from bempy.panel import *
import numpy as np

class TestPanel(unittest.TestCase):
    def test_single_panel_aligned(self):
        x = [1, 0]
        y = [0, 0]
        points = np.vstack([x, y])
        body = Body(points)
        body_panels = BoundVortexPanels(body)
        body_panels.update_strengths()
        xvort, gam = body_panels.vortices
        self.assertEqual(gam, 0)

    def test_single_panel_normal(self):
        x = [0, 0]
        y = [0, 1]
        points = np.vstack([x, y])
        body = Body(points)
        body_panels = BoundVortexPanels(body)
        body_panels.update_strengths()
        xvort, gam = body_panels.vortices
        self.assertEqual(gam, np.pi)
