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
    return s, gam, dgam

num_points = 32
body = TransformedBody(flat_plate(num_points), angle=alpha_deg)
#body = TransformedBody(naca_airfoil('0012',num_points), angle=alpha_deg)
bound = BoundVortices(body)

# unsteady simulation
dt = 0.5
Uinfty = (1,0)
num_steps = 100
stepper = RungeKutta2(dt, Uinfty, bound)

for i in range(num_steps):
    stepper.advance()

print("Using %s timestepper" % stepper.__class__.__name__)
print("Total circulation: %f" % stepper.wake.circulation)
print("After %d steps, last shed vortex has strength: %f" %
      (num_steps, stepper.wake.strengths[-1]))
print("Ratio of last vortex strength to circulation: %f" % 
      (stepper.wake.strengths[-1] / stepper.wake.circulation))
s, gam, dgam = compute_gam(body, stepper.bound.vortices)

pts = body.get_points()
mid = pts[:-1]+0.5*np.diff(pts, axis=0)
xmid = mid[:,0]
rho = 1
qinf = 0.5*rho*(Uinfty[0]**2+Uinfty[1]**2)
# See Kuethe and Chow for the following pressure coefficient formula
cp = -rho*np.linalg.norm(Uinfty)*gam/qinf
print('xmid',xmid)
print('cp',cp)
plt.subplot(211)
plt.plot(xmid,cp,'x')
plt.xlabel('Distance along chord')
plt.ylabel('$C_p$')
plt.title('Computed coefficient of pressure, AoA = %.1f deg' % alpha_deg)

# Exact solution given section 5.4 of Katz/Plotkin, 2nd Ed.
s1 = np.linspace(s[-1],1,100)
dgam_exact = 2 * alpha * np.linalg.norm(Uinfty) * np.sqrt((1-s1) / s1)
plt.subplot(212)
plt.plot(s1, dgam_exact, label="thin airfoil theory")
plt.plot(s, dgam, 's', label="computed, flat plate")
plt.xlabel('Arclength s')
plt.ylabel('$\gamma = d\Gamma/ds$')
plt.title('Comparison with thin airfoil theory, AoA = %.1f deg' % alpha_deg)
plt.grid(True)
plt.legend()
plt.show()
