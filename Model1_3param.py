import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fmin
from scipy.integrate import odeint

# Modèle 1 : Code 2 (3 paramètres)

def dGdIdU (X, t, a, b, Gb):
    Gi = X[0]
    Gs = X[1]
    I = X[2]
    U = X[3]
    dU = - mu * U
    dI = np.maximum(Gs - Gmin, 0) - a * (I - Ib)
    dGi = - U * Gi
    dGs = U * Gi  - b * I * Gs + np.maximum(Gb - Gs, 0)
    return [dGi, dGs, dI, dU]


def xhi2_ogttU (par, t, yd, Gi0):
    b = par[0]
    Gb = par[1]
    U0 = par[2]
    if len(yd) != len(y):
        xhi2 = 1e6
    else :
        xhi2 = sum((y[:, 1] - yd) ** 2) / len(yd)

    if len(par < 0) > 0 :
        xhi2 = xhi2 + 1e5 * len(par < 0)

    if Gb < 70 :
        xhi2 = xhi2 + 1e5
    return xhi2


def compute_absorption (par, Gi0):
    b = par[0]
    Gb = par[1]
    U0 = par[2]
    ya = Gi0 - y[-1, 1]
    return ya


full_ogtt = np.loadtxt("full_ogtt_3param.txt")

tspan = np.linspace(0, 500, 49)
times = full_ogtt[:, 0]

AIC2 = []
params_calib = []
param1 =     [1,    379.804,    826.120,    4.8702E-06,    136.570,    0.003692,
              2,    340.395,    798.721,    4.2955E-06,    121.206,    0.003557,
              3,    1429.144,    703.495,    7.2103E-07,    150,    0.003094,
              4,    522.996,    968.518,    4.6070E-06,    116.613,    0.004411,
              5,    245.002,    767.689,    6.7559E-06,    154.444,    0.003405,
              6,    606.793,    871.109,    4.4321E-06,    85.403,    0.003916,
              7,    368.274,    957.303,    3.6318E-06,    94.607,    0.004354,
              8,    487.455,    981.321,    1.0741E-05,    192.923,    0.004478,
              9,    336.589,    668.115,    3.8786E-06,    171.991,    0.002925,
              10,    340.273,    645.207,    5.3503E-06,    187.893,    0.002817]

[nx, ns] = full_ogtt.shape

for mouse in range(2, (ns - 1)):
    yd = full_ogtt[:, mouse]
    I0 = 10
    Ib = 10
    G0 = yd[0]
    Gmin = G0
    Gb = yd[-1]
    Gi0 = 4500
    a = 1 / 30.0

    mu = 0.0182
    b = 5.15E-006
    U0 = 0.003743

    X0 = [Gi0, G0, I0, U0]
    y = odeint(dGdIdU, X0, tspan, args=(a, b, Gb))

    plt.plot(times, yd, 'go', label="data")
    plt.plot(tspan, y[:, 1], "k-", label="model")
    plt.xlabel('Time(min)')
    plt.ylabel('Glucose(mg/dl)')
    plt.legend()
    plt.show()

    par0 = [0.017422, 5.1632e-06, 0.0038279]
    parOPT = fmin(xhi2_ogttU, par0, args=(tspan, yd, a))
    yabs = compute_absorption(parOPT, Gi0)
    params_calib  = [mouse-1, xhi2_ogttU(parOPT, tspan, yd, Gi0), yabs, parOPT]
    print("params_calib = ", params_calib)
    #sname = 'figs/full_U_{}.pdf'.format(mouse - 1)
    #plt.savefig(sname)