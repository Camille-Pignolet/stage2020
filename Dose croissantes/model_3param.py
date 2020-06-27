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
    ya = Gi0 - y[-1, 0]
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
    elif mouse == 11:
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
plt.plot(times, y1[:, 1], "k-", label="model")
plt.ylabel('')
plt.legend()
plt.title("1g/kg")

plt.subplot(222)
plt.plot(times, data_2g[0], "r.", label="data")
plt.plot(times, data_2g[1], "g.", label="data")
plt.plot(times, data_2g[2], "m.", label="data")
plt.plot(times, y2[:, 1], "k-", label="model")
plt.ylabel('')
plt.legend()
plt.title("2g/kg")

plt.subplot(223)
plt.plot(times, data_3g[0], "r.", label="data")
plt.plot(times, data_3g[1], "g.", label="data")
plt.plot(times, data_3g[2], "m.", label="data")
plt.plot(times, y3[:, 1], "k-", label="model")
plt.xlabel('min')
plt.ylabel('')
plt.legend()
plt.title("3g/kg")

plt.subplot(224)
plt.plot(times, data_4g[0], "r.", label="data")
plt.plot(times, data_4g[1], "g.", label="data")
plt.plot(times, data_4g[2], "m.", label="data")
plt.plot(times, y4[:, 1], "k-", label="model")
plt.xlabel('min')
plt.ylabel('')
plt.legend()
plt.title("4g/kg")
#plt.show()
#plt.close()


#liste des paramètres
#print("params_1g = ", params_doses_1g)
#print("params_2g = ", params_doses_2g)
#print("params_3g = ", params_doses_3g)
#print("params_4g = ", params_doses_4g)

ABS = [47.198, 27.340, 26.354, 16.578]
B = [9.7203E-06, 1.2567E-05, 6.6080E-06, 5.2766E-06]

ET1100 = [54.01, 58.39, 29.19]
ET2100 = [27.34, 34.76, 19.92]
ET3100 = [32.74, 25.02, 21.30]
ET4100 = [19.59, 18.67, 11.47]

SEM1 = np.std(ET1100, axis=0)/np.sqrt(3)
SEM2 = np.std(ET2100, axis=0)/np.sqrt(3)
SEM3 = np.std(ET3100, axis=0)/np.sqrt(3)
SEM4 = np.std(ET4100, axis=0)/np.sqrt(3)

plt.subplot(121)
plt.bar(range(4), ABS, color = ['yellow', 'blue', 'orange', 'pink'], edgecolor = 'black', yerr = [SEM1, SEM2, SEM3, SEM4], ecolor = 'black', capsize = 5)
plt.xticks(range(4), ['1g/kg', '2g/kg', '3g/kg', '4g/kg'])
myText = plt.title("% d'Absorption")
myText.set_fontsize(16)

B1100 = [8.6938E-06, 1.4158E-05, 6.3092E-06]
B2100 = [7.3752E-06, 1.3540E-05, 1.6787E-05]
B3100 = [3.7270E-06, 8.4829E-06, 7.6141E-06]
B4100 = [5.8683E-06, 4.4001E-06, 5.5614E-06]

SEM5 = np.std(B1100, axis=0)/np.sqrt(3)
SEM6 = np.std(B2100, axis=0)/np.sqrt(3)
SEM7 = np.std(B3100, axis=0)/np.sqrt(3)
SEM8 = np.std(B4100, axis=0)/np.sqrt(3)

plt.subplot(122)
plt.bar(range(4), B, color = ['yellow', 'blue', 'orange', 'pink'], edgecolor = 'black', yerr = [SEM5, SEM6, SEM7, SEM8], ecolor = 'black', capsize = 5)
plt.xticks(range(4), ['1g/kg', '2g/kg', '3g/kg', '4g/kg'])
myText1 = plt.title("Sensibilité à l'insuline (min-1)")
myText1.set_fontsize(16)
plt.show()