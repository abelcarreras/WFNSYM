import numpy as np


def rotation_xy(alpha, beta):
    Rx = [[1, 0, 0],
          [0, np.cos(alpha), np.sin(alpha)],
          [0, -np.sin(alpha), np.cos(alpha)]]

    Ry = [[np.cos(beta), 0, -np.sin(beta)],
          [0, 1, 0],
          [np.sin(beta), 0, np.cos(beta)]]

    return np.dot(Rx, Ry)


def target_function_test(alpha, beta):
    c_geom = np.sum(data['coordinates'], axis=0)/len(data['coordinates'])
    ax1 = np.dot(rotation_xy(alpha, beta), [1, 0, 0])
    ax2 = np.dot(rotation_xy(alpha, beta), [0, 0, 1])

    return c_geom, ax1, ax2


def first_approximation(target_function, delta, center=None):

    x = y = np.arange(0.0, np.pi, delta)
    X, Y = np.meshgrid(x, y)

    min = 100
    xmin = 0
    ymin = 0
    for vx, vy in zip(X, Y):
        for xx,yy in zip(vx, vy):
            if center is None:
                val = target_function(xx, yy)
            else:
                val = target_function(xx, yy, center)
            if val < min:
                min = val
                xmin = xx
                ymin = yy

    return xmin, ymin, min


def minimize_axis(target_function, center_i, delta=0.05):
    from scipy.optimize import fmin

    xmin, ymin, val = first_approximation(target_function, delta, center=center_i)

    def minf(x):
        if len(x) > 2:
            center_i = [x[2], x[3], x[4]]
            return target_function(x[0], x[1], center=center_i)
        else:
            return target_function(x[0], x[1])

    if center_i is None:
        return fmin(minf,[xmin, ymin])
    else:
        return fmin(minf,[xmin, ymin, center_i[0], center_i[1], center_i[2]])


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

    print(minimize_axis(target_function, center_i=[0, 0, 0], delta=0.05))
