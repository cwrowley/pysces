import numpy as np

def naca4(code, num_points, zero_thick_te=False, uniform=False):
    """Return coordinates for a NACA 4-digit series airfoil
    """

    # extract parameters from 4-digit code
    code_str = "%04d" % int(code)
    if len(code_str) != 4:
        raise ValueError("Designation is more than 4 digits")
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
    y_thick = 5 * thickness * (np.polyval(coefs[:5], x) + coefs[5] * np.sqrt(x))

    # camber
    front = np.where(x <= p)
    back = np.where(x > p)
    y_camber = np.zeros_like(x)
    if p:
        y_camber[front] = max_camber * x[front] / p**2 * (2 * p - x[front])
        y_camber[back] = max_camber * ((1. - x[back])/(1. - p)**2 *
                                       (1 + x[back] - 2 * p))
    xout = np.hstack([x[-1::-1], x])
    yout = np.hstack([y_camber[-1::-1] - y_thick[-1::-1], y_camber + y_thick])
    return xout, yout

def main():
    import matplotlib.pyplot as plt
    x,y = naca4("0006", 40, zero_thick_te=True)
    x2, y2 = naca4("2406", 40)
    print([(x[0],y[0]),(x[-1],y[-1])])
    plt.plot(x, y, '+-', x2, y2, 'x-')
    plt.axis('equal')
    plt.show()

if __name__ == "__main__":
    main()
