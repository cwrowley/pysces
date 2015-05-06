from pysces import *
import numpy as np
import matplotlib.pyplot as plt

panels = Vortices()
panels.core_radius = 0.01
x = np.linspace(0,0.2,200)
y = np.zeros_like(x)
q = np.array([x,y]).T
xvort = np.array([0,0])
vel = panels.induced_velocity_single(q, xvort, 1)
plt.plot(x,vel[:,0], '-x', label="u")
plt.plot(x,vel[:,1], '-+', label="v")
plt.legend()
plt.show()
