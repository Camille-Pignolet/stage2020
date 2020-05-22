import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fmin
from scipy.integrate import odeint

# Modèle 1 : Code 3 (2 paramètres)

def dGdIdU (X, t, Gna, b, Gb):
    Gi = X[0]
    Gs = X[1]
    I = X[2]
    U = X[3]
    dU = - 0.017 * U
    dI = np.maximum(Gs - Gmin, 0) - (I - Ib) / 30
    dGi = - U * (Gi - Gna)
    dGs = U * (Gi - Gna) - b * I * Gs + np.maximum(Gb - Gs, 0)
    return [dGi, dGs, dI, dU]


def xhi2_ogttU (par, t, yd, Gi0):
    b = par[0]
    U0 = par[1]
    if len(yd) != len(y):
        xhi2 = 1e6
    else :
        xhi2 = sum((y[:, 1] - yd) ** 2) / len(yd)

    if len(par < 0) > 0 :
        xhi2 = xhi2 + 1e5 * len(par < 0)

    if (Gna > Gi0) :
        xhi2 = xhi2 + 1e5

    if Gb < 70 :
        xhi2 = xhi2 + 1e5
    return xhi2


def compute_absorption (par, Gi0, ):
    b = par[0]
    U0 = par[1]
    ya = Gi0 - y[-1, 1]
    return ya


full_ogtt = np.loadtxt("full_ogtt.txt")

tspan = np.linspace(0, 500, 49)
times = full_ogtt[:, 0]

AIC2 = []
params_calib = []
param1 = [1,    384.448,    889.040,     5.15E-006,        0.003743,
          2,    367.177,     833.324,     4.42E-006,        0.003482,
          3,    1480.647,    739.959,     3.32E-006,        0.003055,
          4,    633.270,     988.660,     4.64E-006,        0.004218,
          5,    240.818,     810.507,     6.97E-006,        0.003377,
          6,    720.789,     878.740,     4.42E-006,        0.003694,
          7,    428.898,     966.703,     3.64E-006,        0.004112,
          8,    577.943,     1004.113,    1.07E-005,        0.004293,
          9,    435.298,     778.152,     4.34E-006,        0.003228,
          10,    428.621,     632.897,     5.16E-006,        0.002577]

[nx, ns] = full_ogtt.shape

for mouse in range(2, (ns - 1)):
    yd = full_ogtt[:, mouse]
    I0 = 10
    Ib = 10
    Gna = 0
    G0 = yd[0]
    Gmin = G0
    Gb = yd[-1]
    Gi0 = 4500
    a = 1 / 30.0

    b = 5.15E-006
    U0 = 0.003743

    X0 = [Gi0, G0, I0, U0]
    y = odeint(dGdIdU, X0, tspan, args=(Gna, b, Gb))

    plt.plot(times, yd, 'r.', label="data")
    plt.plot(tspan, y[:, 1], "b-", label="model")
    plt.xlabel('Time(min)')
    plt.ylabel('Glucose(mg/dl)')
    plt.legend()
    plt.show()

    par0 = [10, 0.002577]
    parOPT = fmin(xhi2_ogttU, par0, args=(tspan, yd, a))
    yabs = compute_absorption(parOPT, Gi0)
    params_calib  = [mouse-1, xhi2_ogttU(parOPT, tspan, yd, Gi0), yabs, parOPT]
    AIC2 = [49 * np.log(xhi2_ogttU(parOPT, tspan, yd, a) / 49) + 6]
    print(params_calib)
    print(AIC2)
    #sname = 'figs/full_U_{}.pdf'.format(mouse - 1)
    #plt.savefig(sname)