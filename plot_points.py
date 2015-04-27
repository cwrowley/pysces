from bempy import *
import matplotlib.pyplot as plt

airfoil = Airfoil("0012", 6)   # NACA 2412 airfoil with 20 points per side

body_panels = BoundVortexPanels(airfoil)

xv, yv, gv = body_panels.vortices
xc, yc = body_panels.collocation_pts
x, y = airfoil.get_points()
plt.plot(x, y, 'k-+')
plt.plot(xv, yv, 'ro')
plt.plot(xc, yc, 'bx')
plt.axis('equal')
plt.grid(True)
plt.show()
