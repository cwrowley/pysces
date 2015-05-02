import numpy as np
from pysces import *

# x = [1,0]
# y = [0,0]
# x = [0,0,0]
# y = [0,1,2]
x = [2,1,0]
y = [0,0,0]
points = np.array([x, y]).T
body = Body(points)
body_panels = BoundVortexPanels(body)
body_panels.update_strengths()
gam = body_panels.vortices.strengths
xc = body_panels.collocation_pts
Uinfty = np.array((1,0))
vel = body_panels.induced_velocity(xc)
vel_dot_n = np.sum(vel * body_panels.normals, 0)
print("gam:" + str(gam))
print("vel: " + str(vel))
print("vel . n:" + str(vel_dot_n))
