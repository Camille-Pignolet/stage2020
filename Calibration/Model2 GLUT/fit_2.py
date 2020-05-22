import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fmin
from scipy.integrate import odeint

# Modèle 2 : Code 2

def dGdIdUdGlut (X, t, a, b, mu, Gb, c, s, Lambda):
    Gi = X[0]
    Gs = X[1]
    I = X[2]
    U = X[3]
    Glut = X[4]
    dU = - mu * U
    dI = np.maximum(Gs - Gmin, 0) - a * (I - Ib)
    dGi = - U * (1+Glut) * Gi
    dGs = U * (1+Glut) * Gi - b * I * Gs + np.maximum(Gb - Gs, 0)
    dGlut = U * Lambda * ((Gs - G0) / (Gs + s)) - c * (I - Ib) * Glut
    if Glut < 0:
        dGlut = 0
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

    if Lambda < 100 :
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


full_ogtt = np.loadtxt("full_ogtt2.txt")

tspan = np.linspace(0, 500, 49)
times = full_ogtt[:, 0]

params_calib = []
param1 = [1,	172.993,	526.857,	2.9627E-006,	106.606,	0.00131175,	0.00267145,	10684.6900,
        2,	121.684,	518.229,	2.8190E-006,	150.000,	0.00116446,	0.00058231,	4814.7827,
        3,	1050.545,	567.226,	2.9734E-006,	186.328,	0.00120315,	0.00003694,	831.0899,
        4,	388.480,	789.260,	3.9875E-006,	142.182,	0.00264846,	0.00007881,	167.8535,
        5,	255.033,	679.409,	6.1487E-006,	158.737,	0.00271818,	0.00023436,	100.0000,
        6,	419.444,	637.632,	3.3240E-006,	94.616,	0.00178914,	0.00027162,	1002.2284,
        7,	414.079,	822.405,	3.1612E-006,	121.198,	0.00316581,	0.00014543,	100.0000,
        8,	387.769,	578.918,	6.6709E-006,	190.171,	0.00176722,	0.00368942,	6427.4831,
        9,	296.419,	571.731,	3.1895E-006,	160.956,	0.00214925,	0.00320845,	1441.2956,
        10,	306.812,	439.610,	3.7191E-006,	185.633,	0.00134154,	0.00057795,	1574.2669]


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
    s = 200
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

    par0 = [4e-06, yd[-1], 0.001, 0.0005, 2000]
    parOPT = fmin(xhi2_ogttU, par0, args=(tspan, yd, Gi0, s))
    yabs = compute_absorption(parOPT, Gi0)
    params_calib  = [mouse-1, xhi2_ogttU(parOPT, tspan, yd, Gi0, s), yabs, parOPT]
    print(params_calib)
    #sname = 'figs/full_U_{}.pdf'.format(mouse - 1)
    #plt.savefig(sname)

  