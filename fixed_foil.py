from bempy import *
import matplotlib.pyplot as plt

airfoil = Airfoil("2412", 20)   # NACA 2412 airfoil with 20 points per side
airfoil = TransformedBody(airfoil, displacement=(-0.25, 0))
airfoil = TransformedBody(airfoil, angle=10) # rotate by 10 degrees about 1/4 chord

body_panels = BoundVortexPanels(airfoil)
wake_panels = FreeVortexPanels()

num_steps = 10
Uinfty = (1,0)
dt = 0.1
for i in range(num_steps):
    time = i * dt
    time_advance(body_panels, wake_panels, Uinfty, dt)
    lift, drag = compute_forces(body_panels, wake_panels)
    print("Time %.1f: Lift = %.3f, Drag = %.3f" % (time, lift, drag))

xv, yv, gv = wake_panels.vortices
x, y = airfoil.get_points()
plt.plot(x, y, 'k-+')
plt.plot(xv, yv, 'ro')
plt.axis('equal')
plt.grid(True)
plt.show()
