import numpy as np
import matplotlib.pyplot as plt # - это граф библиотека.
import math
import scipy.integrate as spi
import random
import numpy as np
import matplotlib.pyplot as plt
import math
from tkinter import Tk, Label
import pandas as pd

def Monte_Curve(te, te1, te2, N, N1, N2, percentbadf, NeededAmount, dt): # это я создаю распределение. Не суть важно
 Stretch=te2/te1*10*percentbadf/100
 if percentbadf<=1:
     Stretch=50
 te1=te1*Stretch
 te2=te2*Stretch  
 s1=te1/((te2-te1)/2)**0.5
 s2=te1/((te2-te1)/2)**0.5
 
 s1=(te2-te1)*0.05/(math.log(2))**0.5
 s1=1/(2*math.pi)**0.5/0.5
 s2=s1

 g11 =lambda x: (np.exp(-1*(x-te1)**2/(2*s1**2)))*(1/((math.pi*2)**0.5*s1))
 from scipy import integrate
 Res11=integrate.quad(g11, 0, np.inf)[0]



 maxx=1.1*te2

 g3= lambda x: (((np.exp(-1*(x-te1)**2/(2*s1**2)))*(1/((math.pi*2)**0.5*s1))*N1/N)/Res11+((np.exp(-1*(x-te2)**2/(2*s2**2)))*(1/((math.pi*2)**0.5*s2))*N2/N))
 from scipy import integrate
 Res=integrate.quad(g3, -np.inf, np.inf)[0]
 #print("интеграл g3=",Res)
 g4= lambda x: (((np.exp(-1*(x-te1)**2/(2*s1**2)))*(1/((math.pi*2)**0.5*s1))*N1/N)/Res11+((np.exp(-1*(x-te2)**2/(2*s2**2)))*(1/((math.pi*2)**0.5*s2))*(N2/N)))/Res
 Resg4=integrate.quad(g4, 0, np.inf)[0]
 #print("интеграл g4=",Resg4)
 g5= lambda x: (((np.exp(-1*(x-te1)**2/(2*s1**2)))*(1/((math.pi*2)**0.5*s1))*N1/N)/Res11+((np.exp(-1*(x-te2)**2/(2*s2**2)))*(1/((math.pi*2)**0.5*s2))*(N2/N)))/Resg4
 #Resg5=integrate.quad(g5, 0, maxx)[0]
 #print("интеграл финальной функц. распределения g5=",Resg5)
 #g5x= lambda x: (((np.exp(-1*(x-te1)**2/(2*s1**2)))*(1/((math.pi*2)**0.5*s1))*N1/N)+((np.exp(-1*(x-te2)**2/(2*s2**2)))*(1/((math.pi*2)**0.5*s2))*(N2/N)))/Resg4*x
 #MidWeight=integrate.quad(g5x, 0, maxx)[0]/Resg5
 #print("Среднее Взвешанное =",MidWeight/Stretch, "при пике te1 =",te1/Stretch, "при пике te2 =",te2/Stretch)
 #x = np.linspace(0, 1.1*te2 ,1000000)
 #average = sum (x) / len (x)
 #print("Среднее =",average/Stretch)
 
 """
 plt.plot(x, g3(x))
 plt.show()
 

 plt.plot(x, g4(x))
 plt.show()

 plt.plot(x, g5(x))
 plt.show()
 """
 Maxg=g5(te1)
 Ming=0
 minx=0
 Teapr=[]
 Gteapr=[]
 i=0
 while i <= NeededAmount:
  Condition = False
  while Condition !=True:
   X=random.uniform(minx, maxx)
   Y=random.uniform(Ming, Maxg)
   if Y <= g5(X):
    Condition = True
    teApproved=X
    gteApproved=Y
    Teapr.append(teApproved)
    Gteapr.append(gteApproved)
    i=i+1
 dtnew = 0.1*min(Teapr)
 if dt > dtnew:
     dt2=dtnew
 else:
     dt2=dt
 Teapr=np.array(Teapr)
 Tbarier=(te2-te1)/2
 N1gone= Teapr[Teapr<=Tbarier].shape[0]
 N2gone=NeededAmount-N1gone
 x=np.linspace(0, te2*1.1, 10000)
 plt.plot(x, g5(x))
 plt.show()
 Picke=[g5(te1), g5(te2)]
 pickedots=[te1, te2]
 plt.plot(x, g5(x))
 plt.plot(Teapr, Gteapr, '.')
 plt.plot(pickedots, Picke, 'o')
 plt.show()

 return Teapr, dt2 , N1gone, N2gone



V = 2500*10**(-9)*10**9
Ve = 100*10**(-18)*10**9

Np = int(10**5)*6
L = 4*10**(-6)*10**3 #длина поры
d = 3.5*10**(-6)*10**3 #диаметр поры
a = 1/1.41374728*10**3/32 # относительный параметр, отвечающий за скорость течения жидкости через одну пору
b = a/100*8.55 #относительный параметр, отвечающий за скорость продвижения эритроцита через мембрану (относительно скорости воды "а")

te = (L + 8*Ve/(math.pi*d**2))/b
percentnonf=0#Процент нефильтруемых
percentbadf=5

Gemotocrid=1
N = int(round(Gemotocrid/100*V/Ve))
N1= int(round(Gemotocrid/100*V/Ve))*(1-percentbadf/100)
N2= int(round(Gemotocrid/100*V/Ve))*(percentbadf/100)


te1 =te #-исправить, это должно быть мат ожидание, от нуля до 1
te2 = te1*100

NeededAmount=1200
npe=N
dt=te1/10
Teapr, dt2, N1gone, N2gone = Monte_Curve(te, te1, te2, npe, N1, N2, percentbadf, NeededAmount, dt)
print(N1gone, N2gone)
print((N2gone/(N1gone+ N2gone)*100))