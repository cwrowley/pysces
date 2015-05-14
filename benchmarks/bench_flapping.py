import numpy as np
from pysces import *
from timeit import default_timer as timer

start = timer()
airfoil = naca_airfoil("2412", 20)   # NACA 2412 airfoil with 20 points per side
airfoil = TransformedBody(airfoil, displacement=(-0.25, 0))
freq = 0.3 * 2 * np.pi
airfoil = Pitching(airfoil, 10, freq, phase=90)
airfoil = Heaving(airfoil, (0,0.2), freq, phase=0)
bound = BoundVortices(airfoil)

num_steps = 400
Uinfty = (1,0)
dt = 0.01
Vortices.core_radius = dt

stepper = RungeKutta2(dt, Uinfty, bound)

print("Taking %d steps" % num_steps)
for i in range(1,num_steps):
    stepper.advance()

elapsed = timer() - start
print("Time: %f sec" % elapsed)
