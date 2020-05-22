import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fmin
from scipy.integrate import odeint

# Modèle 2 : Code 1

def dGdIdUdGlut (X, t, a, b, mu, Gb, c, s, Lambda):
    Gi = X[0]
    Gs = X[1]
    I = X[2]
    U = X[3]
    Glut = X[4]
    dU = - mu * U * (1+Glut)
    dI = np.maximum(Gs - Gmin, 0) - a * (I - Ib)
    dGi = - U * Gi - Glut * Gi
    dGs = U * Gi + Glut * Gi - b * I * Gs + np.maximum(Gb - Gs, 0)
    dGlut = U * Lambda * (Gs / (Gs * s)) - c * (I - Ib) * Glut
    return [dGi, dGs, dI, dU, dGlut]


def xhi2_ogttU (par, t, yd, Gi0, s):
    b = par[0]
    Gb = par[1]
    U0 = par[2]
    c = par[3]
    Lambda = par[4]
    if len(yd) != len(y):
        xhi2 = 1e6
    else :
        xhi2 = sum((y[:, 1] - yd) ** 2) / len(yd)

    if len(par < 0) > 0 :
        xhi2 = xhi2 + 1e5 * len(par < 0)

    if Gb < 70 :
        xhi2 = xhi2 + 1e5

    if c < 1e-7 :
        xhi2 = xhi2 + 1e5
    return xhi2


def compute_absorption (par, Gi0, ):
    b = par[0]
    Gb = par[1]
    U0 = par[2]
    c = par[3]
    Lambda = par[4]
    ya = Gi0 - y[-1, 1]
    return ya


full_ogtt = np.loadtxt("full_ogtt1.txt")

tspan = np.linspace(0, 500, 49)
times = full_ogtt[:, 0]

params_calib = []
param1 = [1,    206.313,    635.545,    3.4431E-006,    108.549,    0.00185873,    0.00053949,    0.4121,
          2,    158.424,    568.900,    2.8345E-006,    124.429,    0.00088153,    0.00013061,    0.4817,
          3,    1138.879,    853.292,    3.5516E-006,    100.158,    0.00351750,    0.00002408,    0.01,
          4,    353.035,    1057.902,    4.8295E-006,    122.433,    0.00316122,    0.00002129,    0.0386,
          5,    330.837,    515.945,    4.6727E-006,    169.570,    0.00056324,    0.00009287,    0.6345,
          6,    410.953,    689.322,    3.3524E-006,    70.301,    0.00140949,    0.00008435,    0.2426,
          7,    367.717,    957.336,    3.6331E-006,    77.519,    0.00435382,    0.00012580,    0.0000,
          8,    488.889,    948.743,    1.0359E-005,    193.133,    0.00411400,    0.00027544,    0.0236,
          9,    345.359,    448.468,    2.5456E-006,    174.781,    0.00080662,    0.00012419,    0.4677,
          10,    322.007,    481.369,    3.9966E-006,    188.692,    0.00132903,    0.00017984,    0.2392]


[nx, ns] = full_ogtt.shape

for mouse in range(2, (ns - 1)):
    yd = full_ogtt[:, mouse]
    I0 = 10
    Ib = 10
    Glut0 = 0
    Gna = 0
    G0 = yd[0]
    Gmin = G0
    Gi0 = 4500
    a = 1 / 30.0
    s = 1.1
    mu = 0.0182

    b = 3.5516E-006
    Gb = 100
    U0 = 0.00351750
    c = 0.00002408
    Lambda = 0.01

    X0 = [Gi0, G0, I0, U0, Glut0]
    y = odeint(dGdIdUdGlut, X0, tspan, args=(a, b, mu, Gb, c, s, Lambda))

    plt.subplot(231)
    plt.plot(times, yd, 'r.', label="data")
    plt.plot(tspan, y[:, 1], "b-", label="model")
    plt.xlabel('Time(min)')
    plt.ylabel('Glucose(mg/dl)')
    plt.legend()
    plt.subplot(232)
    plt.plot(tspan, y[:, 0], "b-", label="model")
    plt.xlabel('Time(min)')
    plt.ylabel('Gi')
    plt.legend()
    plt.subplot(233)
    plt.plot(tspan, y[:, 1], "b-", label="model")
    plt.xlabel('Time(min)')
    plt.ylabel('Gs')
    plt.legend()
    plt.subplot(234)
    plt.plot(tspan, y[:, 2], "b-", label="model")
    plt.xlabel('Time(min)')
    plt.ylabel('I')
    plt.legend()
    plt.subplot(235)
    plt.plot(tspan, y[:, 3], "b-", label="model")
    plt.xlabel('Time(min)')
    plt.ylabel('U')
    plt.legend()
    plt.subplot(236)
    plt.plot(tspan, y[:, 4], "b-", label="model")
    plt.xlabel('Time(min)')
    plt.ylabel('Glut')
    plt.legend()

    plt.show()

    par0 = [3e-06, yd[-1], 0.002, 0.0001, 0.4]
    parOPT = fmin(xhi2_ogttU, par0, args=(tspan, yd, Gi0, s))
    yabs = compute_absorption(parOPT, Gi0)
    params_calib  = [mouse-1, xhi2_ogttU(parOPT, tspan, yd, Gi0, s), yabs, parOPT]
    print(params_calib)
    #sname = 'figs/full_U_{}.pdf'.format(mouse - 1)
    #plt.savefig(sname)



