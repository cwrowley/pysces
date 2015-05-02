from pysces import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

airfoil = naca_airfoil("0006", 20) # NACA 0012 airfoil with 20 points per side
airfoil = TransformedBody(airfoil, displacement=(-0.25, 0))
freq = 0.3 * 2*np.pi
airfoil = Pitching(airfoil, 20, freq, phase=90)
airfoil = Heaving(airfoil, (0,0.2), freq, phase=0)

Uinfty = (1,0)
dt = 0.05
Vortices.core_radius = dt
flow = Simulation(airfoil, Uinfty, dt, BoundVortices)

fig, ax = plt.subplots()
ax.axis('equal')
ax.axis([-1,5,-2,2])
ax.grid(True)
q = airfoil.get_points()
line, = ax.plot(q[:,0], q[:,1], '-k')
pts, = ax.plot(0,0,'ro')

def gen_points():
    flow.initialize()
    num_steps = 200
    dt = 0.02
    for i in range(num_steps):
        flow.advance()
        yield airfoil.get_points(), flow.wake.positions

def redraw(data):
    q, xvort = data
    line.set_data(q[:,0], q[:,1])
    pts.set_data(xvort[:,0], xvort[:,1])

movie = animation.FuncAnimation(fig, redraw, gen_points, interval=50, repeat_delay=0)
plt.show()
