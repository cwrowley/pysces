from pysces import *
import matplotlib.pyplot as plt

airfoil = naca_airfoil("2412", 20)   # NACA 2412 airfoil with 20 points per side
airfoil = TransformedBody(airfoil, displacement=(-0.25, 0))
airfoil = TransformedBody(airfoil, angle=10) # rotate by 10 degrees about 1/4 chord

num_steps = 100
Uinfty = (1,0)
dt = 0.01
Vortices.core_radius = dt

flow = ExplicitEuler(dt, Uinfty, airfoil, BoundVortices)
flow2 = RungeKutta2(dt, Uinfty, airfoil, BoundVortices)

for i in range(1,num_steps):
    flow.advance()
    flow2.advance()
    # lift, drag = compute_forces(flow.body_panels, flow.wake_panels)
    # print("Time %.1f: Lift = %.3f, Drag = %.3f" % (flow.time, lift, drag))

vort = flow.wake.positions
vort2 = flow2.wake.positions
q = airfoil.get_points()
plt.plot(q[:,0], q[:,1], 'k-')
plt.plot(vort[:,0], vort[:,1], 'ro')
plt.plot(vort2[:,0], vort2[:,1], 'go')
plt.axis('equal')
plt.grid(True)
plt.show()
