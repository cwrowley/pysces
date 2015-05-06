from pysces import *
import matplotlib.pyplot as plt
import numpy as np

v1 = (-2,0)
v2 = (1,0)
g1 = 2 * np.pi
g2 = 4 * np.pi
vort = Vortices([v1,v2],[g1,g2])
omega = 1./3

dt = 0.2
Uinfty = (0,0)
euler = ExplicitEuler(dt, Uinfty, wake=vort)
rk2 = RungeKutta2(dt, Uinfty, wake=vort)
rk4 = RungeKutta4(dt, Uinfty, wake=vort)

num_steps = 10
t = dt * np.arange(num_steps)

def exact(t):
    v = np.array([np.cos(omega * t), np.sin(omega * t)])
    return np.array([-2 * v, v])

for i in range(1,num_steps):
    euler.advance()
    rk2.advance()
    rk4.advance()
    q1 = euler.wake.positions
    q2 = rk2.wake.positions
    q3 = rk4.wake.positions
    v = exact(t[i])
    err = [np.sum((q - v) * (q - v)) for q in (q1, q2, q3)]
    print(err)
    plt.plot(q1[:,0], q1[:,1], 'r+')
    plt.plot(q2[:,0], q2[:,1], 'b+')
    plt.plot(q3[:,0], q3[:,1], 'g+')
    plt.plot(v[:,0], v[:,1], 'x')

plt.axis('equal')
plt.grid(True)
plt.show()
