from __future__ import division
import numpy as np
from pysces import *
import matplotlib.pyplot as plt

alpha_deg = 2 # degrees
alpha = alpha_deg * np.pi / 180

def compute_gam(body, vort):
    q = body.get_points()
    gam = vort.strengths
    ds = -np.linalg.norm(np.diff(q, axis=0), axis=1)
    dgam = gam / ds
    xvort = vort.positions
    s = np.sqrt(xvort[:,0]**2 + xvort[:,1]**2)
    return s, dgam

num_points = 32
angle = 2 # degrees
body = flat_plate(num_points)
body = TransformedBody(body, angle=angle)
bound = BoundVortices(body)

# unsteady simulation
dt = 0.5
Uinfty = (1,0)
num_steps = 100
# stepper = ExplicitEuler(dt, Uinfty, bound)
stepper = RungeKutta2(dt, Uinfty, bound)

for i in range(num_steps):
    stepper.advance()

vort = stepper.wake.positions
print("Using %s timestepper" % stepper.__class__.__name__)
print("Total circulation: %f" % stepper.wake.circulation)
print("After %d steps, last shed vortex: %f" %
      (num_steps, stepper.wake.strengths[-1]))
print("Ratio: %f" % (stepper.wake.strengths[-1] / stepper.wake.circulation))
s, dgam = compute_gam(body, stepper.bound.vortices)

s1 = np.linspace(s[-1],1,100)
dgam_exact = 2 * alpha * np.sqrt((1-s1) / s1)

plt.plot(s1, dgam_exact, label="thin airfoil theory")
plt.plot(s, dgam, 'x', label="computed, flat plate")
plt.xlabel('Arclength s')
plt.ylabel(r'$\gamma = d\Gamma/ds$')
plt.title('Comparison with thin airfoil theory, AoA = %.1f deg' % alpha_deg)
plt.grid(True)
plt.legend()
plt.show()
