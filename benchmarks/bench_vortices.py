import numpy as np
from pysces import Vortices
from timeit import default_timer as timer

n = 8192
Vortices.core_radius = 0.01

pos = np.array(np.random.rand(2*n), dtype=np.float32).reshape((n,2))
gamma = np.array(np.random.rand(n), dtype=np.float32)
print("Computing induced velocity for %d vortices" % n)

start = timer()
vort = Vortices(pos, gamma)
vel = vort.induced_velocity()
elapsed = timer() - start

print("Time: %f sec" % elapsed)
