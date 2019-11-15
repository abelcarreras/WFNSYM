import numpy as np
from scipy.optimize import fmin


def rotation_xy(alpha, beta):
    Rx = [[1, 0, 0],
          [0, np.cos(alpha), np.sin(alpha)],
          [0, -np.sin(alpha), np.cos(alpha)]]

    Ry = [[np.cos(beta), 0, -np.sin(beta)],
          [0, 1, 0],
          [np.sin(beta), 0, np.cos(beta)]]

    return np.dot(Rx, Ry)


# Rotation respect to axis a
def rotation_axis(a, angle):
    cos_1 = 1 - np.cos(angle)
    R = [[a[0]*a[0]*cos_1 + np.cos(angle),      a[0]*a[1]*cos_1 - a[2]*np.sin(angle), a[0]*a[2]*cos_1 + a[1]*np.sin(angle)],
         [a[1]*a[0]*cos_1 + a[2]*np.sin(angle), a[1]*a[1]*cos_1 + np.cos(angle),      a[1]*a[2]*cos_1 - a[0]*np.sin(angle)],
         [a[2]*a[0]*cos_1 - a[1]*np.sin(angle), a[2]*a[1]*cos_1 + a[0]*np.sin(angle), a[2]*a[2]*cos_1 + np.cos(angle)]]
    return np.array(R)


def target_function_test(alpha, beta):
    c_geom = np.sum(data['coordinates'], axis=0)/len(data['coordinates'])
    ax1 = np.dot(rotation_xy(alpha, beta), [1, 0, 0])
    ax2 = np.dot(rotation_xy(alpha, beta), [0, 0, 1])

    return c_geom, ax1, ax2


def first_approximation(target_function, delta, min_gamma=False):

    x = y = np.arange(0.0, np.pi, delta)
    X, Y = np.meshgrid(x, y)
    z = np.arange(0.0, np.pi, delta)

    min = None
    xmin = 0
    ymin = 0
    zmin = 0
    for vx, vy in zip(X, Y):
        for xx,yy in zip(vx, vy):
            if min_gamma:
                for zz in z:
                    val = target_function(xx, yy, gamma=zz)
                    if min is None or val < min:
                        min = val
                        xmin = xx
                        ymin = yy
                        zmin = zz
            else:
                val = target_function(xx, yy, gamma=0.0)
                if  min is None or val < min:
                    min = val
                    xmin = xx
                    ymin = yy
                    zmin = 0.0
                    # print(xx, yy, val)

    return xmin, ymin, zmin, min


def first_approximation_gamma_only(target_function, vaxis, delta):
    z = np.arange(0.0, np.pi, delta)

    min = None
    zmin = 0
    for zz in z:
        val = target_function(zz, vaxis)
        if min is None or val < min:
            min = val
            zmin = zz
            # print(xx, yy, zz, val)

    return zmin, min


def minimize_axis(target_function, data, delta=0.05):

    # define optimization functions
    def minf(x):
        return target_function(x[0], x[1])

    def minf_g(x):  # with gamma
        return target_function(x[0], x[1], gamma=x[2])

    # Check if second axis is needed
    minimize_ax2 = True if data['igroup'] == 8 else False

    xmin, ymin, zmin, _ = first_approximation(target_function, delta, min_gamma=minimize_ax2)
    if minimize_ax2:
        alpha, beta, gamma = fmin(minf_g, [xmin, ymin, zmin], disp=False)
        return alpha, beta, gamma
    else:
        alpha, beta = fmin(minf, [xmin, ymin], disp=False)
        return alpha, beta, 0


def minimize_axis2(target_function, axis, delta=0.05):

    # define optimization functions
    def minf(x):
        return target_function(x[0], axis)

    zmin, _ = first_approximation_gamma_only(target_function, axis, delta)

    gamma = fmin(minf, [zmin], disp=False)
    return gamma[0]


if __name__ == "__main__":
    from wfnsympy.file_io import get_data_from_file_fchk
    from wfnsympy import WfnSympy

    data = get_data_from_file_fchk('../unittest/pirrol.in.fchk')
    c_geom = np.sum(data['coordinates'], axis=0)/len(data['coordinates'])

    def target_function(alpha, beta):

        ax1 = np.dot(rotation_xy(alpha, beta), [1, 0, 0])
        # ax2 = np.dot(rotation_xy(alpha, beta), [0, 0, 1])

        pirrol = WfnSympy(coordinates=data['coordinates'],
                          symbols=data['symbols'],
                          basis=data['basis'],
                          center=c_geom,
                          VAxis=ax1,
                          # VAxis2=ax2,
                          alpha_mo_coeff=data['mo_coefficients']['alpha'],
                          charge=0,
                          multiplicity=1,
                          group='C5')
        res = np.sum(np.sqrt(pirrol.csm_coef))
        return res

    print(minimize_axis(target_function, data={'igroup': 2}, delta=0.05))
