# from typing import Any, Union

import numpy as np

def tspline(x,y, n, xout, yout, iout, sigma, yp, temp):
    '''
    :param x: double
    :param y: double
    :param n: integer
    :param xout: double
    :param iout: int
    :param sigma: float
    :param yp: double
    :param tmp: double
    :return:
    '''
    nm1 = n-2
    np1 = n
    delx1 = x[1] - x[0]
    dx1 = (y[1] - y[0]) / delx1
    delx2 = x[2] - x[1]
    delx12 = x[2] - x[0]
    c1 = -(delx12 + delx1) / delx12 / delx1
    c2 = delx12 / delx1 / delx2
    c3 = -delx1 / delx12 / delx2
    slpp1 = c1 * y[0] + c2 * y[1] + c3 * y[2]
    deln = x[n - 1] - x[nm1]
    delnm1 = x[nm1] - x[n - 3]
    delnn = x[n - 1] - x[n - 3]
    c1 = (delnn + deln) / delnn / deln
    c2 = -delnn / deln / delnm1
    c3 = deln / delnn / delnm1
    slppn = c3 * y[n - 3] + c2 * y[nm1] + c1 * y[n - 1]
    sigmap = np.fabs(sigma) * (float(n - 2)) / (x[n - 1] - x[0])
    dels = sigmap * delx1
    exps = np.exp(dels)
    sinhs = 0.5 * (exps - 1. / exps)
    sinhin = 1. / (delx1 * sinhs)
    diag1 = sinhin * (dels * 0.5 * (exps + 1. / exps) - sinhs)
    diagin = 1. / diag1
    yp[0] = diagin * (dx1 - slpp1)
    spdiag = sinhin * (sinhs - dels)
    temp[0] = diagin * spdiag
    dx2 = 0
    for i in range(0, nm1):
        delx2 = x[i+1] - x[i]
        dx2 = (y[i+1]-y[i]) / delx2
        dels = sigmap * delx2
        exps = np.exp(dels)
        sinhs = 0.5 * (exps-1. / exps)
        sinhin = 1. / (delx2 * sinhs)
        diag2 = sinhin * (dels * (0.5 * (exps+1. / exps))-sinhs)
        diagin = 1. / (diag1+diag2-spdiag * temp[i-1])
        yp[i] = diagin * (dx2-dx1-spdiag * yp[i-1])
        spdiag = sinhin * (sinhs-dels)
        temp[i] = diagin * spdiag
        dx1 = dx2
        diag1 = diag2
    diagin = 1. / (diag1 - spdiag * temp[nm1])
    yp[n - 1] = diagin * (slppn - dx2 - spdiag * yp[nm1])
    for i in range(1, n):
        ibak = np1 - i - 1
        yp[ibak] = yp[ibak] - temp[ibak] * yp[ibak + 1]
    a = x[0]
    b = x[1]
    nj = 1
    for i in range(0,iout):
        if (xout[i] > b):
            while (xout[i] > b):
                a = b
                nj = nj+1
                b = x[nj]

        del1 = xout[i] - a
        del2 = b - xout[i]
        dels = b - a
        exps1 = np.exp(sigmap * del1)
        sinhd1 = 0.5 * (exps1 - 1. / exps1)
        exps = np.exp(sigmap * del2)
        sinhd2 = 0.5 * (exps - 1. / exps)
        exps = exps * exps1
        sinhs = 0.5 * (exps - 1. / exps)
        yout[i] = (yp[nj] * sinhd1 + yp[nj - 1] * sinhd2) / \
                  sinhs + ((y[nj] - yp[nj]) * del1 + (y[nj - 1] - yp[nj - 1]) * del2) / dels