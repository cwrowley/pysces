from bempy import *
import matplotlib.pyplot as plt

airfoil = Airfoil("0012", 20)      # NACA 0012 airfoil with 20 points per side
airfoil = TransformedBody(airfoil, displacement=(-0.25, 0))
airfoil = Pitching(airfoil, 10, 2*np.pi, phase=90)
airfoil = Heaving(airfoil, (0,0.2), 2*np.pi, phase=0)

body_panels = SourceDoubletPanels(airfoil)
wake_panels = FreeVortexParticles()

num_steps = 10
Uinfty = (1,0)
dt = 0.1
for i in range(num_steps):
    time = i * dt
    airfoil.time = time
    body_panels.update_positions()
    time_advance(body_panels, wake_panels, Uinfty, dt)
    # some sort of output for body position and wake vortex positions/strengths
    lift, drag = compute_forces(body_panels, wake_panels)
    print("Time %.1f: Lift = %.3f, Drag = %.3f" % (time, lift, drag))

xv, yv, gv = wake_panels.vortices
x, y = airfoil.get_points()
plt.plot(x, y, 'k-+')
plt.plot(xv, yv, 'ro')
plt.axis('equal')
plt.grid(True)
plt.show()
