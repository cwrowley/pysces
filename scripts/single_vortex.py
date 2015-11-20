from pysces import *
import matplotlib.pyplot as plt
import numpy as np

v1 = (0,0)
v2 = (1,0)
s1 = 2 * np.pi
s2 = 0
vort = Vortices([v1,v2],[s1,s2])

dt = 0.1
Uinfty = (0,0)
euler = ExplicitEuler(dt, Uinfty, wake=vort)
rk2 = RungeKutta2(dt, Uinfty, wake=vort)
rk4 = RungeKutta4(dt, Uinfty, wake=vort)

num_steps = 200
for i in range(1,num_steps):
    euler.advance()
    rk2.advance()
    rk4.advance()
    q1 = euler.wake.positions
    q2 = rk2.wake.positions
    q3 = rk4.wake.positions
    plt.plot(q1[:,0], q1[:,1], 'ro')
    plt.plot(q2[:,0], q2[:,1], 'bo')
    plt.plot(q3[:,0], q3[:,1], 'go')

plt.axis('equal')
plt.grid(True)
plt.show()
