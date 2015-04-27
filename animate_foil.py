from bempy import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation

airfoil = Airfoil("0012", 20)      # NACA 0012 airfoil with 20 points per side
airfoil = TransformedBody(airfoil, displacement=(-0.25, 0))
airfoil = Pitching(airfoil, 10, 2*np.pi, phase=90)
airfoil = Heaving(airfoil, (0,0.2), 2*np.pi, phase=0)
airfoil = TransformedBody(airfoil, displacement=(3,4))

fig, ax = plt.subplots(1)
ax.axis('equal')
ax.grid(True)
q = airfoil.get_points()
line, = ax.plot(q[0], q[1], '-k+')

def gen_points():
    num_steps = 50
    dt = 0.02
    for i in range(num_steps):
        airfoil.time = i * dt
        yield airfoil.get_points()

def redraw(q):
    line.set_data(q[0], q[1])

movie = animation.FuncAnimation(fig, redraw, gen_points, interval=50, repeat_delay=0)
plt.show()
