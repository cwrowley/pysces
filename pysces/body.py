import numpy as np
from .motion import RigidMotion

__all__ = ['Body', 'TransformedBody', 'Pitching', 'Heaving',
           'cylinder', 'flat_plate', 'naca_airfoil', 'joukowski_foil',
           'van_de_vooren_foil', 'karman_trefftz_foil']

class Body(object):
    """Base class for representing bodies
    """
    def __init__(self, points):
        """Create a body with nodes at the given points

        Parameters
        ----------
        points : 2d array, shape (n,2)
            Array of points defining the boundary of the body
            For a closed body, the boundary curve should be positively oriented
            (counter-clockwise around outside of body), starting from trailing
            edge
        """
        self._time = 0
        self._points = np.array(points, dtype="float64")

    @property
    def time(self):
        """The time used to specify the body's motion"""
        return self._time

    @time.setter
    def time(self, value):
        self._time = value

    def get_points(self, body_frame=False):
        return self._points

    def get_body(self):
        """Return the Body object in the body-fixed frame"""
        return self

    def get_motion(self):
        """Return the transformation from the body-fixed to inertial frame"""
        return None

def cylinder(radius, num_points):
    """Return a circular Body with the given radius and number of points"""
    th = np.linspace(0, 2 * np.pi, num_points)
    points = radius * np.array([np.cos(th), np.sin(th)]).T
    return Body(points)

def flat_plate(num_points):
    """Return a flat plate Body with the given number of points.
    
    In body coordinates the plate runs from (0,0) to (1,0)."""
    x = np.linspace(1, 0, num_points) # from 1 to 0 so trailing edge at index 0
    y = np.zeros_like(x)
    return Body(np.array([x, y]).T)

def joukowski_foil(xcenter=-.1, ycenter=.1, a=1, numpoints=32):
    """Return a Joukowski foil Body.
    
    The foil has its trailing edge at (2a,0).  The foil has a total of
    numpoints along the boundary.  Refer to chapter 4 of [1]_ for details.

    Parameters
    ----------
    xcenter, ycenter : float
        (xcenter,ycenter) is the center of the Joukowski preimage circle.
        xcenter should be negative and small; its magnitude determines
        the bluffness of the foil.  ycenter should be small; it
        determines the magnitude of the camber (positive gives upward
        camber, and negative gives downward camber).

    a : float
        radius of the Joukowski preimage circle

    numpoints : int
        number of points along the boundary

    References
    ----------
    .. [1] Acheson, D. J., "Elementary Fluid Dynamics", Oxford, 1990.
    """

    t = np.linspace(0,2*np.pi,numpoints)
    r = np.sqrt((a-xcenter)**2+ycenter**2)
    chi = xcenter + r*np.cos(t)
    eta = ycenter + r*np.sin(t)
    mag2 = chi*chi + eta*eta
    x = chi*(1+a**2/mag2)
    y = eta*(1-a**2/mag2)
    return Body(np.array([x,y]).T)

def karman_trefftz_foil(xcenter=-.1, ycenter=0, a=.1, angle_deg=10, numpoints=32):
    """Return a Karman-Trefftz foil Body.
    
    The Karman-Trefftz foil is a modified version of the Joukowski
    foil but with a nonzero interior angle --- rather than a cusp ---  at the 
    trailing edge.  Refer to [1]_ for details.
    
    Parameters
    ----------
    xcenter, ycenter, a : float.
        The same as in joukowski_foil().

    angle_deg : float
        The interior angle, in degrees, at the trailing edge.

    numpoints : int
        Number of points along the boundary

    See Also
    --------
    joukowski_foil()

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Joukowsky_transform
    """

    angle_rad = angle_deg*np.pi/180
    n = 2-angle_rad/np.pi
    t = np.linspace(0,2*np.pi,numpoints)
    ctr = xcenter + 1j*ycenter
    r = np.linalg.norm(ctr-a)
    zeta = ctr+r*np.exp(1j*t)
    mag2 = np.linalg.norm(zeta)
    z = n*((1+1/zeta)**n+(1-1/zeta)**n)/((1+1/zeta)**n-(1-1/zeta)**n)
    x = [w.real for w in z]
    y = [w.imag for w in z]
    return Body(np.array([x,y]).T)

def van_de_vooren_foil(semichord=1.0, thickness=0.15, angle_deg=5,
numpoints=32):
    """Return a van de Vooren foil Body.

    Refer to section 6.6 of [1]_

    Parameters
    ----------
    semichord : float
        half the chord c, so c=2*semichord

    thickness : float
        vertical thickness as a fraction (0 < thickness < 1) of the semichord

    angle_deg : float
        interior angle, in degrees, at the trailing edge

    numpoints : int
        number of points along the boundary

    References
    ----------
    .. [1] Katz, Joseph and Plotkin, Allen, "Low-Speed Aerodynamics", 2nd Ed.,
       Cambridge University Press, 2001.
    """

    k = 2-(angle_deg*np.pi/180)
    a = 2*semichord*((1+thickness)**(k-1))*2**(-k)
    t = np.linspace(0,2*np.pi,numpoints)
    num = (a*(np.cos(t)-1)+1j*a*np.sin(t))**k
    den = (a*(np.cos(t)-thickness)+1j*a*np.sin(t))**(k-1)
    z = (num/den)+semichord
    x = [w.real for w in z]
    y = [w.imag for w in z]
    return Body(np.array([x,y]).T)

def naca_airfoil(code, num_points, zero_thick_te=False, uniform=False):
    """Return a NACA 4-digit series airfoil"""
    # extract parameters from 4-digit code
    code_str = "%04d" % int(code)
    if len(code_str) != 4:
        raise ValueError("NACA designation is more than 4 digits")
    max_camber = 0.01 * int(code_str[0])
    p = 0.1 * int(code_str[1])  # location of max camber
    thickness = 0.01 * int(code_str[2:])
    if uniform:
        x = np.linspace(0, 1, num_points)
    else:
        # closer spacing near leading edge
        theta = np.linspace(0, 0.5 * np.pi, num_points)
        x = 1 - np.cos(theta)

    # thickness
    coefs = [-0.1015, 0.2843, -0.3516, -0.1260, 0, 0.2969]
    if zero_thick_te:
        coefs[0] = -0.1036
    y_thick = 5 * thickness * (np.polyval(coefs[:5], x) +
                               coefs[5] * np.sqrt(x))

    # camber
    front = np.where(x <= p)
    back = np.where(x > p)
    y_camber = np.zeros_like(x)
    if p:
        y_camber[front] = max_camber * x[front] / p**2 * (2 * p - x[front])
        y_camber[back] = max_camber * ((1. - x[back])/(1. - p)**2 *
                                       (1 + x[back] - 2 * p))
    x = np.hstack([x[-1:0:-1], x])
    y = np.hstack([y_camber[-1:0:-1] + y_thick[-1:0:-1],
                   y_camber - y_thick])
    return Body(np.array([x, y]).T)


class TransformedBody(object):
    """Base class for rigid (Euclidean) transformations of existing bodies
    """
    def __init__(self, body, angle=0, displacement=(0,0)):
        """angles are clockwise, in degrees"""
        self._parent = body
        self._body = body.get_body()
        self._motion = RigidMotion(-angle * np.pi / 180, displacement)

    def get_body(self):
        return self._body

    def get_motion(self):
        self._update()
        return self._motion.compose(self._parent.get_motion())

    def set_motion(self, value):
        self._motion = value

    @property
    def time(self):
        return self._body.time

    @time.setter
    def time(self, value):
        self._body.time = value

    def _update(self):
        # update body motion: subclasses override this
        pass

    def get_points(self, body_frame=False):
        q = self._body.get_points()
        if body_frame:
            return q
        return self.get_motion().map_position(q)


class Pitching(TransformedBody):
    """Sinusoidal pitching for an existing body
    """
    def __init__(self, body, amplitude, frequency, phase=0.):
        """amplitude and phase given in degrees"""
        super(Pitching, self).__init__(body)
        self._amplitude = amplitude * np.pi / 180
        self._frequency = frequency
        self._phase = phase * np.pi / 180

    def _update(self):
        theta = self._frequency * self.time + self._phase
        alpha = self._amplitude * np.sin(theta)
        alphadot = self._amplitude * self._frequency * np.cos(theta)
        self.set_motion(RigidMotion(-alpha, (0,0), -alphadot, (0,0)))


class Heaving(TransformedBody):
    """Sinusoidal heaving for an existing body
    """
    def __init__(self, body, displacement, frequency, phase=0.):
        super(Heaving, self).__init__(body)
        self._displacement = np.array(displacement, dtype="float64")
        self._frequency = frequency
        self._phase = phase * np.pi / 180

    def _update(self):
        theta = self._frequency * self.time + self._phase
        x = self._displacement * np.sin(theta)
        xdot = self._displacement * self._frequency * np.cos(theta)
        self.set_motion(RigidMotion(0, x, 0, xdot))
