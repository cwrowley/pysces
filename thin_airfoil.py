from __future__ import division
import numpy as np
from pysces import *
import matplotlib.pyplot as plt

alpha_deg = 2 # degrees
alpha = alpha_deg * np.pi / 180

def compute_gam(body):
    body = TransformedBody(body, angle=alpha_deg, displacement=(1,2))
    panels = BoundVortices(body)
    panels.update_strengths()
    xvort, gam = panels.vortices.positions, panels.vortices.strengths
    q = body.get_points()
    ds = np.linalg.norm(np.diff(q, axis=0), axis=1)
    # import pdb; pdb.set_trace()
    dgam = gam / ds
    s = np.sqrt(xvort[:,0]**2 + xvort[:,1]**2)
    return s, dgam

num_points = 32
plate = flat_plate(num_points)
airfoil = naca_airfoil("0001", num_points)

s_plate, dgam_plate = compute_gam(plate)

s_airfoil, dgam_airfoil = compute_gam(airfoil)
# sum up vortices on top and bottom of airfoil
half = dgam_airfoil.shape[0] // 2
dgam_airfoil = dgam_airfoil[:half] + dgam_airfoil[-1:half-1:-1]
s_airfoil = s_airfoil[:half]

# exact distribution from thin airfoil theory: see Kuethe and Chow p143?
s = np.linspace(s_plate[-1],1,100)
dgam_exact = 2 * alpha * np.sqrt((1-s) / s)

plt.plot(s, dgam_exact, label="thin airfoil theory")
plt.plot(s_plate, dgam_plate, 'x', label="computed, flat plate")
plt.plot(s_airfoil, dgam_airfoil, '+', label="computed, NACA 0001")
plt.xlabel('Arclength s')
plt.ylabel(r'$\gamma = d\Gamma/ds$')
plt.ylim([0,1])
plt.title('Comparison with thin airfoil theory, AoA = %.1f deg' % alpha_deg)
plt.grid(True)
plt.legend()
plt.show()
