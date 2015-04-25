from bempy import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation

airfoil = Airfoil("0012", 20)      # NACA 0012 airfoil with 20 points per side
airfoil = Translation(airfoil, displacement=(-0.25, 0))
airfoil = Translation(airfoil, angle=0)
airfoil = Pitching(airfoil, 10, 2*np.pi, phase=90)
airfoil = Heaving(airfoil, (0,0.2), 2*np.pi, phase=0)

fig, ax = plt.subplots(1)
ax.axis('equal')
ax.grid(True)
x, y = airfoil.get_points()
line, = ax.plot(x, y, '-k+')

def gen_points():
    num_steps = 50
    dt = 0.02
    for i in range(num_steps):
        airfoil.time = i * dt
        x, y = airfoil.get_points()
        yield x, y

def redraw(data):
    x, y = data
    line.set_data(x, y)

movie = animation.FuncAnimation(fig, redraw, gen_points, interval=50, repeat_delay=0)
plt.show()
