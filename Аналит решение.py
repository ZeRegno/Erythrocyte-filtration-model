
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.integrate as spi
import random
import numpy as np
import matplotlib.pyplot as plt
import math
from tkinter import Tk, Label
import pandas as pd

V = 2500*10**(-9)*10**9
Ve = 100*10**(-18)*10**9

Np = int(int(10**5)*6)
Np2 = int(int(10**5)*6)*2
L = 4*10**(-6)*10**3 #длина поры
d = 3.5*10**(-6)*10**3 #диаметр поры
a = 1/1.41374728*10**3/32 # относительный параметр, отвечающий за скорость течения жидкости через одну пору
b = a/100*8.55 #относительный параметр, отвечающий за скорость продвижения эритроцита через мембрану (относительно скорости воды "а")

percentnonf=0#Процент нефильтруемых

percentbadf=10
Percentdouble=0
Gemotocrid=5

te = (L + 8*Ve/(math.pi*d**2))/b
 # время за которое эритроцит проходит через пору
te=round(te, 6)
te1=te
te2=te*100
dt = 10**(-1)*te

v0 =a*d**2/4*math.pi
N0 = int(round(Gemotocrid/100*V/Ve))

te=te

Vвых1=1/te1*10**(-3)

Vвых12=1/te1*10**(-3)*2
Vзакр1=Gemotocrid/100-Gemotocrid/100*percentbadf/100

Vвых2=1/te2*10**(-1)
Vзакр2=Gemotocrid/100*percentbadf/100

dV1= lambda t: -(-v0*Np+(v0-Vвых1*Ve)*Np*Vзакр1/(Vзакр1+Vвых1)*(-np.exp(-Vзакр1*t-Vвых1*t)+1))/v0/Np
dV12= lambda t: -(-v0*Np+(v0-Vвых12*Ve)*Np*Vзакр1/(Vзакр1+Vвых12)*(-np.exp(-Vзакр1*t-Vвых12*t)+1))/v0/Np

dV12Np= lambda t: -(-v0*Np2+(v0-Vвых1*Ve)*Np2*Vзакр1/(Vзакр1+Vвых1)*(-np.exp(-Vзакр1*t-Vвых1*t)+1))/v0/Np2

V= lambda t: -(-v0*Np*t+(v0-Vвых1*Ve)*Np*Vзакр1/(Vзакр1+Vвых1)*(np.exp(-Vзакр1*t-Vвых1*t)/(Vзакр1+Vвых1)+t))/v0/Np

V0= lambda t: t*0

dV2= lambda t: (-v0*Np+((Np*Vзакр1/(Vзакр1+Vвых1)-(np.exp(-Vзакр1*t-Vвых1*t)+1))*(v0-Vзакр1*Vвых1))+((Np*Vзакр2/(Vзакр2+Vвых2)-(np.exp(-Vзакр2*t-Vвых2*t)+1))*(v0-Vзакр2*Vвых2)))/(-v0*Np)

t = np.linspace(0, 100 ,1000000)

plt.plot(t, dV1(t))
plt.plot(t, dV12Np(t))
plt.plot(t, V0(t))

plt.show()

"""
plt.plot(t, dV1(t))
plt.plot(t, V0(t))
plt.show()
"""