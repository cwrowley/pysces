from bempy import *
import matplotlib.pyplot as plt

airfoil = Airfoil("8412", 6)   # NACA 2412 airfoil with 20 points per side
airfoil = TransformedBody(airfoil, angle=10)

# cylinder = Circle(1.0, 21)

body_panels = BoundVortexPanels(airfoil)

xv, yv, gv = body_panels.vortices
xc, yc = body_panels.collocation_pts
xn, yn = body_panels.normals
x, y = airfoil.get_points(body_frame=True)
plt.plot(x, y, 'k-+')
plt.plot(xv, yv, 'ro', label="vortices")
plt.plot(xc, yc, 'bx', label="collocation pts")
base = np.vstack([xc,yc])
plt.quiver(xc, yc, xn, yn)
plt.legend()
plt.axis('equal')
plt.grid(True)
plt.show()
