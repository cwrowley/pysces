from __future__ import division
import numpy as np
from bempy import *
import matplotlib.pyplot as plt

body = naca_airfoil("0012", 16)
body = flat_plate(16)
body = TransformedBody(body, angle=10)
panels = BoundVortexPanels(body)
Uinfty = (1,0)
dt = 1
panels.update_strengths_unsteady(None, Uinfty, dt)

_, gam = panels.vortices
x_shed, gam_shed = panels.get_newly_shed()
gam_sum = np.sum(gam)
print("gam_shed = %f, bound circ = %f, total = %f" % (gam_shed, gam_sum, gam_shed+gam_sum))

# plot
q = body.get_points()
plt.plot(q[0], q[1], 'k-')
plt.plot(x_shed[0], x_shed[1], 'bo')
plt.grid(True)
plt.axis('equal')
plt.show()
