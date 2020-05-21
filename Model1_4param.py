import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fmin
from scipy.integrate import odeint

# Modèle 1 : Code 1 (4 paramètres)

def dGdIdU (X, t, Gna, a, b, mu, Gb) :
    Gi = X[0]
    Gs = X[1]
    I = X[2]
    U = X[3]
    Ib = 10
    Gmin = G0
    dU = - mu * U
    dI = np.maximum(Gs-Gmin, 0) - a * (I-Ib)
    dGi = - U * (Gi - Gna)
    dGs = - U * (Gi-Gna) - b * I * Gs + np.maximum(Gb-Gs, 0)
    return [dGi, dGs, dI, dU]


def xhi2_ogttU(par, t, yd, a) :
    mu = par[0]
    b = par[1]
    Gb = par[2]
    U0 = par[3]
    if len(yd) != len(y) :
        xhi2 = 1e6
    else :
        xhi2 = sum((y[:, 1] - yd)**2) / len(yd)

    if len(par < 0) > 0 :
        xhi2 = xhi2 + 1e5 * len(par < 0)

    if (Gna > Gi0) :
        xhi2 = xhi2 + 1e5

    if Gb < 70 :
        xhi2 = xhi2 + 1e5
    return xhi2


def compute_absorption(par, Gi0 ) :
    mu = par[0]
    b = par[1]
    Gb = par[2]
    U0 = par[3]
    ya = Gi0 - y[-1, 1]
    return ya


full_ogtt = np.loadtxt("full_ogtt_4param.txt")

tspan = np.linspace(0, 500, 49)
times = full_ogtt[:, 0]

AIC1 = []
params = []
param1 = [379.04,       853.92,     0.017465,   4.9835e-06,       137.49,    0.0036757,
          339.64,       824.72,       0.0175,   4.3965e-06,        99.56,    0.0035435,
          1321.8,        564.7,     0.024777,   2.7098e-06,       165.03,    0.0033225,
          519.31,       918.49,     0.019479,   4.4412e-06,       115.43,    0.0044472,
          236.33,       873.24,     0.015603,   7.2761e-06,       155.49,    0.0033676,
          599.9,       949.34,    0.016358,   4.7078e-06,       93.099,    0.0038771,
          360.34,       862.61,     0.020832,   3.3475e-06,       73.706,    0.0044335,
          449.76,       1308.6,     0.013055,   1.2142e-05,       191.53,    0.0044924,
          282.64,       1134.2,     0.010094,   5.3048e-06,       165.69,    0.0029505,
          357.59,       1076.3,    0.0098727,   6.7358e-06,       188.25,    0.0027183]

[nx, ns] = full_ogtt.shape

for mouse in range(2, (ns - 1)) :
    yd = full_ogtt[:, mouse]
    I0 = 10
    Ib = 10
    Gna = 0
    G0 = yd[0]
    Gmin = G0
    Gi0 = 4500
    a = 1 / 30.0

    mu = 0.017
    b = 4.9835e-06
    Gb = 151
    U0 = 0.0036757

    X0 = [Gi0, G0, I0, U0]
    y = odeint(dGdIdU, X0, tspan, args=(Gna, a, b, mu, Gb,))

    plt.plot(times, yd, 'r.', label="data")
    plt.plot(tspan, y[:, 1], "b-", label="model")
    plt.xlabel('Time(min)')
    plt.ylabel('Glucose(mg/dl)')
    plt.legend()
    plt.show()

    par0 = [0.017422, 5.1632e-06, 70.306, 0.0038279]
    parOPT = fmin(xhi2_ogttU, par0, args=(times, yd, a))
    yabs = compute_absorption(parOPT, Gi0)
    params = [xhi2_ogttU(parOPT, tspan, yd, Gi0), yabs, parOPT]
    AIC1 = [49 * np.log(xhi2_ogttU(parOPT, times, yd, a) / 49) + 8]
    print(params)
    print(AIC1)
    # sname = 'figs/full_U_{}.pdf'.format(mouse - 1)
    # plt.savefig(sname)