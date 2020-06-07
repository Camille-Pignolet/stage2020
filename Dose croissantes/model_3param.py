import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fmin
from scipy.integrate import odeint
import csv


# Modèle 1 : équations, xhi2, compute absorption
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
    y = odeint(dGdIdU, X0, tspan, args=(a, b, Gb))
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
    y = odeint(dGdIdU, X0, tspan, args=(a, b, Gb))
    ya = Gi0 - y[-1, 1]
    return ya


# Import data dose croissantes
with open('OGTT Doses croissantes.csv', newline='') as csvfile:
    dose_croissantes = csv.reader(csvfile, delimiter=';')

    times = []
    Dose_data = []
    next(dose_croissantes)
    next(dose_croissantes)
    for row in dose_croissantes:
        times.append(int(row[0]))
        for column in range(1, len(row)-1, 3):
            Dose_data.append((int(row[column]), int(row[column+1]), int(row[column+2])))

    data_1g = ([], [], [])
    data_2g = ([], [], [])
    data_3g = ([], [], [])
    data_4g = ([], [], [])
    dose_data_index = 0
    for x in range(len(times)):
        data_1g[0].append(Dose_data[dose_data_index][0])
        data_1g[1].append(Dose_data[dose_data_index][1])
        data_1g[2].append(Dose_data[dose_data_index][2])

        data_2g[0].append(Dose_data[dose_data_index+1][0])
        data_2g[1].append(Dose_data[dose_data_index+1][1])
        data_2g[2].append(Dose_data[dose_data_index+1][2])

        data_3g[0].append(Dose_data[dose_data_index + 2][0])
        data_3g[1].append(Dose_data[dose_data_index + 2][1])
        data_3g[2].append(Dose_data[dose_data_index + 2][2])

        data_4g[0].append(Dose_data[dose_data_index + 3][0])
        data_4g[1].append(Dose_data[dose_data_index + 3][1])
        data_4g[2].append(Dose_data[dose_data_index + 3][2])

        dose_data_index += 4


# scrip to fit model
tspan = np.linspace(0, 180, 49)
params_doses_1g = []
params_doses_2g = []
params_doses_3g = []
params_doses_4g = []

param1 = [1, 701.882, 593.822, 8.6792E-06, 194.781, 0.01413,
      2, 794.338, 642.349, 1.4157E-05, 220.166, 0.01596,
      3, 224.128, 321.091, 6.3092E-06, 78.857, 0.00628,
      4, 456.551, 601.538, 7.3752E-06, 109.040, 0.00581,
      5, 1010.966, 766.977, 1.3554E-05, 200.659, 0.00780,
      6, 1101.780, 440.344, 1.6829E-05, 250.056, 0.00407,
      7, 386.188, 1078.973, 3.7225E-06, 169.888, 0.00721,
      8, 1113.936, 825.796, 8.4829E-06, 78.675, 0.00524,
      9, 668.195, 702.844, 7.6141E-06, 211.322, 0.00436,
      10, 316.423, 861.955, 5.8695E-06, 183.933, 0.00397,
      11, 769.925, 821.697, 4.4009E-06, 240.273, 0.00376,
      12, 809.367, 504.726, 5.5614E-06, 252.410, 0.00222]

I0 = 10
Ib = 10
a = 1 / 30.0
mu = 0.0182
b = 5.15E-006
U0 = 0.003743

for mouse in range(0, 2):
    if mouse == 0 :
        yd = data_1g[0]
    elif mouse == 1 :
        yd = data_1g[1]
    elif mouse == 2 :
        yd = data_1g[2]
    G0 = yd[0]
    Gmin = G0
    Gb = yd[-1]
    Gi0 = 1100
    X0 = [Gi0, G0, I0, U0]
    y1 = odeint(dGdIdU, X0, tspan, args=(a, b, Gb))
    par0 = [5E-06, yd[-1], 0.004]
    parOPT = fmin(xhi2_ogttU, par0, args=(tspan, yd, a))
    yabs = compute_absorption(parOPT, Gi0)
    params_doses_1g = [mouse - 1, xhi2_ogttU(parOPT, tspan, yd, Gi0), yabs, parOPT]


for mouse in range(3, 5):
    if mouse == 3:
        yd = data_2g[0]
    elif mouse == 4:
        yd = data_2g[1]
    elif mouse == 5:
        yd = data_2g[2]
    G0 = yd[0]
    Gmin = G0
    Gb = yd[-1]
    Gi0 = 2200
    X0 = [Gi0, G0, I0, U0]
    y2 = odeint(dGdIdU, X0, tspan, args=(a, b, Gb))
    par0 = [5E-06, yd[-1], 0.004]
    parOPT = fmin(xhi2_ogttU, par0, args=(tspan, yd, a))
    yabs = compute_absorption(parOPT, Gi0)
    params_doses_2g = [mouse - 1, xhi2_ogttU(parOPT, tspan, yd, Gi0), yabs, parOPT]


for mouse in range(6, 8):
    if mouse == 6:
        yd = data_3g[0]
    elif mouse == 7:
        yd = data_3g[1]
    elif mouse == 8:
        yd = data_3g[2]
    G0 = yd[0]
    Gmin = G0
    Gb = yd[-1]
    Gi0 = 3300
    X0 = [Gi0, G0, I0, U0]
    y3 = odeint(dGdIdU, X0, tspan, args=(a, b, Gb))
    par0 = [5E-06, yd[-1], 0.004]
    parOPT = fmin(xhi2_ogttU, par0, args=(tspan, yd, a))
    yabs = compute_absorption(parOPT, Gi0)
    params_doses_3g = [mouse - 1, xhi2_ogttU(parOPT, tspan, yd, Gi0), yabs, parOPT]


for mouse in range(9, 11):
    if mouse == 9:
        yd = data_4g[0]
    elif mouse == 10:
        yd = data_4g[1]
    elif mouse == 21:
        yd = data_4g[2]
    G0 = yd[0]
    Gmin = G0
    Gb = yd[-1]
    Gi0 = 4400
    X0 = [Gi0, G0, I0, U0]
    y4 = odeint(dGdIdU, X0, tspan, args=(a, b, Gb))
    par0 = [5E-06, yd[-1], 0.004]
    parOPT = fmin(xhi2_ogttU, par0, args=(tspan, yd, a))
    yabs = compute_absorption(parOPT, Gi0)
    params_doses_4g = [mouse - 1, xhi2_ogttU(parOPT, tspan, yd, Gi0), yabs, parOPT]


# plot des fits
plt.subplot(221)
plt.plot(times, data_1g[0], "r.", label="data")
plt.plot(times, data_1g[1], "g.", label="data")
plt.plot(times, data_1g[2], "m.", label="data")
plt.plot(times, y1[:, 1], "-", label="model")
plt.ylabel('')
plt.legend()
plt.title("1g/kg")

plt.subplot(222)
plt.plot(times, data_2g[0], "r.", label="data")
plt.plot(times, data_2g[1], "g.", label="data")
plt.plot(times, data_2g[2], "m.", label="data")
plt.plot(times, y2[:, 1], "-", label="model")
plt.ylabel('')
plt.legend()
plt.title("2g/kg")

plt.subplot(223)
plt.plot(times, data_3g[0], "r.", label="data")
plt.plot(times, data_3g[1], "g.", label="data")
plt.plot(times, data_3g[2], "m.", label="data")
plt.plot(times, y3[:, 1], "-", label="model")
plt.xlabel('min')
plt.ylabel('')
plt.legend()
plt.title("3g/kg")

plt.subplot(224)
plt.plot(times, data_4g[0], "r.", label="data")
plt.plot(times, data_4g[1], "g.", label="data")
plt.plot(times, data_4g[2], "m.", label="data")
plt.plot(times, y4[:, 1], "-", label="model")
plt.xlabel('min')
plt.ylabel('')
plt.legend()
plt.title("4g/kg")
plt.show()
plt.close()


#liste des paramètres
print("params_1g = ", params_doses_1g)
print("params_2g = ", params_doses_2g)
print("params_3g = ", params_doses_3g)
print("params_4g = ", params_doses_4g)

XHI2 = (params_doses_1g[1]+params_doses_2g[1]+params_doses_3g[1]+params_doses_4g[1])/4
YABS = (params_doses_1g[2]+params_doses_2g[2]+params_doses_3g[2]+params_doses_4g[2])/4
B = (params_doses_1g[3][0]+params_doses_2g[3][0]+params_doses_3g[3][0]+params_doses_4g[3][0])/4
MU = (params_doses_1g[3][1]+params_doses_2g[3][1]+params_doses_3g[3][1]+params_doses_4g[3][1])/4
UO = (params_doses_1g[3][2]+params_doses_2g[3][2]+params_doses_3g[3][2]+params_doses_4g[3][2])/4

plt.subplot(121)
plt.bar(range(2), [XHI2, YABS])
plt.xticks(range(2), ['XHI2', 'yabs'])

plt.subplot(122)
plt.bar(range(3), [B, MU, UO])
plt.xticks(range(3), ['b', 'MU', 'U0'])
plt.show()