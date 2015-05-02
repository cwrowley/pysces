from pysces import *
import matplotlib.pyplot as plt
import numpy as np

airfoil = naca_airfoil("0012", 20)      # NACA 0012 airfoil with 20 points per side
airfoil = TransformedBody(airfoil, displacement=(-0.25, 0))
freq = 0.7
airfoil = Pitching(airfoil, 10, freq, phase=90)
airfoil = Heaving(airfoil, (0,0.2), freq, phase=0)

num_steps = 100
Uinfty = (1,0)
dt = 0.1

flow = Simulation(airfoil, Uinfty, dt, BoundVortexPanels)

for i in range(1,num_steps):
    flow.advance()

vort = flow.wake.positions
q = airfoil.get_points()
plt.plot(q[:,0], q[:,1], 'k-')
plt.plot(vort[:,0], vort[:,1], 'ro')
plt.axis('equal')
plt.grid(True)
plt.show()
