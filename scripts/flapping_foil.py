from pysces import *
import matplotlib.pyplot as plt
import numpy as np

airfoil = naca_airfoil("0012", 50)      # NACA 0012 airfoil with 20 points per side
airfoil = TransformedBody(airfoil, displacement=(-0.25, 0))
freq = 0.3 * 2 * np.pi
airfoil = Pitching(airfoil, 10, freq, phase=90)
airfoil = Heaving(airfoil, (0,0.2), freq, phase=0)

num_steps = 400
Uinfty = (1,0)
dt = 0.01
Vortices.core_radius = dt

flow = ExplicitEuler(dt, Uinfty, airfoil, BoundVortices)

for i in range(1,num_steps):
    flow.advance()

vort = flow.wake.positions
gam = flow.wake.strengths
q = airfoil.get_points()
plt.plot(q[:,0], q[:,1], 'k-')
maxval = dt
plt.scatter(vort[:,0], vort[:,1], c=gam,
            cmap='bwr', vmin=-maxval, vmax=maxval, edgecolors='none')
plt.axis('equal')
plt.grid(True)
plt.show()
