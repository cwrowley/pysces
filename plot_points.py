from bempy import *
import matplotlib.pyplot as plt

airfoil = Airfoil("8412", 6)   # NACA 2412 airfoil with 20 points per side
airfoil = TransformedBody(airfoil, angle=10)

# cylinder = Circle(1.0, 21)

body_panels = BoundVortexPanels(airfoil)

vort, gv = body_panels.vortices
coll = body_panels.collocation_pts
norm = body_panels.normals
foil = airfoil.get_points(body_frame=True)
plt.plot(foil[0], foil[1], 'k-+')
plt.plot(vort[0], vort[1], 'ro', label="vortices")
plt.plot(coll[0], coll[1], 'bx', label="collocation pts")
plt.quiver(coll[0], coll[1], norm[0], norm[1])
plt.legend()
plt.axis('equal')
plt.grid(True)
plt.show()
