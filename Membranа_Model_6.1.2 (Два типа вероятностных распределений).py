import numpy as np
import matplotlib.pyplot as plt
import math
from tkinter import Tk, Label
import pandas as pd

import scipy.integrate as spi
import random




label  =  Label(text = 'Симуляция Мембраны')
label.pack()
window  =  Tk()
def GetSImDATA(V, Np, Gemotocridin, percent, percentbadf, Percentdouble,  Tsim, dt, a , te, te2):
    Gemotocrid=Gemotocridin
    name="Гематокрит_"+ str(round(Gemotocrid, 3))+"Процент нефильтруемых_"+str(round(percent, 3))+"Процент плохофильтр._"+str(round(percentbadf, 3))
    name2="Гематокрит_"+ str(round(Gemotocrid, 3))+"Процент нефильтруемых_"+str(round(percent, 3))+"Процент плохофильтр._"+str(round(percentbadf, 3))
    VEXIT1, t1, dV1, AEM1, Gemotocrit1, AET1, tdV1, dVc1, tGemo1, Gemotocritc1, VEXIT2, t2, dV2, Gemotocrit2, PotokBufer2, PotokEr2, AET2, AEM2, GoodErinMembr2, BadErinMembr2= Sim(V, Np, Gemotocrid, percent, percentbadf, Percentdouble,  Tsim, dt, a , te, te2, name)
    graphs(VEXIT1, t1, dV1, AEM1, Gemotocrit1, AET1, tdV1, dVc1, tGemo1, Gemotocritc1, name2) 
    save(t2, VEXIT2, dV2, Gemotocrit2, PotokBufer2, PotokEr2, AET2, AEM2, GoodErinMembr2, BadErinMembr2, name)

def timesplit(A, t, name): 
 th=1
 n = 1*th
 A1 = []
 T1 = []
 for i in range(int(len(t))):
    T1.append(t[i])
    A1.append(A[i])
    if t[i]//th == n:
        n = n + th
        NAME = str(name) +  "от " + str((n - 2*th)*th) + " сек. до " + str((round((n - 1*th)))*th) + " сек."
        plt.plot(T1, A1,'.', color="red") 
        plt.grid()
        plt.xlabel("Время в с",fontsize = 20)
        plt.title(str(NAME), fontsize = 10)
        plt.figure(figsize = (10000, 5000))
        plt.show()
        A1 = []
        T1 = []
    if max(t) == t[i]:
        NAME = str(name) + " от " + str(n - 1*th) + " сек. до " + str(round(max(t), 3)) + " сек."
        plt.plot(T1, A1,'.', color="red") 
        plt.grid()
        plt.xlabel("Время в с",fontsize = 20)
        plt.title(str(NAME), fontsize = 10)
        plt.figure(figsize = (10000, 5000))
        plt.show()
def cutzero(A, t):
    A1=[]
    T1=[]
    for i in range((len(A))):
     if A[i] != 0:
      if A[i] != max(A):
       if A[i] != min(A):  
         T1.append(t[i])
         A1.append(A[i]) 
    return A1, T1  


def black_bars(expectation, bins, dt):
        poisson = np.random.poisson(expectation, bins)
        bounded_poisson = np.clip(poisson, a_min=0, a_max=dt/2)
        return  bounded_poisson
def black_bars2(expectation, bins):
        poisson = np.random.poisson(expectation, bins)
        bounded_poisson = np.clip(poisson, a_min=0, a_max=1)
        Amount=bounded_poisson[bounded_poisson ==1].shape[0]
        return Amount  
def Monte_Curve(te, te1, te2, N, N1, N2, percentbadf, NeededAmount):
 if percentbadf !=0:
  Stretch=te2/te1*10*percentbadf/100

 if percentbadf<=1:
     Stretch=50
 else: Stretch=10
 te1=te1*Stretch
 te2=te2*Stretch  
 
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
 Maxg=g5(te1)
 Ming=0
 minx=0
 Teapr=[]
 Gteapr=[]
 i=0
 while i <= NeededAmount-1:
  Condition = False
  while Condition !=True:
   X=random.uniform(minx, maxx)
   Y=random.uniform(Ming, Maxg)
   if Y <= g5(X):
    Condition = True
    teApproved=X
    gteApproved=Y
    Teapr.append(teApproved/Stretch)
    Gteapr.append(gteApproved)
    i=i+1
 Teapr=np.array(Teapr)
 Tbarier=((te2-te1)/Stretch)/2
 N1gone= Teapr[Teapr <=Tbarier].shape[0]
 N2gone=NeededAmount-N1gone

 return Teapr, N1gone, N2gone

def Sim(V, Np, Gemotocrid, percent, percentbadf, Percentdouble,  Tsim, dt, a , te, te2, name):
 dVwater =a*d**2/4*math.pi*dt*(Np)/dt
 Ve = 100*10**(-18) *10**9#обьем эритроцита
 Ne = int(round(Gemotocrid/100*V/Ve))
 print("Гематокрит", Gemotocrid)
 print("время прохождения эритроцитом поры", te)
 print("время прохождения плохофильтр. эритроцитом поры", te2)
 Percentdouble=Percentdouble/100

 
 # так тоже читаемо
 VEXIT, dV, t, AEM, AET ,Gemotocrit, PotokBufer, PotokEr, BadErinMembr, GoodErinMembr= (np.zeros(Tsim) for _ in range(10))
 nfe=int(percent/100*Ne) # число нефильтруемых
 print("отношение нефильтруемых к порам", nfe/Np)
 nfeNum = np.linspace(1, int(Ne), int(nfe))
 for i in range(len((nfeNum))):
    nfeNum[i]=round(nfeNum[i])

 npe=Ne # еще не прошедшие через мембрану клетки
 N1=Ne*(1-percentbadf/100)
 N2=Ne*(percentbadf/100)
 p = 0 # количество уже прошедших эритроцитов через мембрану
 ae = 0 # количество активных эритроцитов в мембране 

 blockperc=0
 St  = np.zeros(Np) #массив отвеающий за состояние каждой поры (открыта - 0, занята - любое отличное от нуля число) 
 i = 1 # начальный шаг
 Voverpor=0 
 Vwater = 0 # обьем воды прошедшей через мембрану
 Vwaterall = V - Ne*Ve # обьем всей воды участвующей в расчетах, в не зависимости от состояния
 GemoVEXIT=Gemotocrid
 Gemostart=Gemotocrid
 """
 я высчитывал пройдет или непройдет эритроцит через мембрану при помощи таких этих порций
 суть в следующем, если данный обьем порции кратен целому числу обьемов воды,
 содержащих только 1 эритроцит (это число я считал ниже исходя из равномерной распределенности эритроцитов в жидкости)
 то считается сколько таких целых обьемов уместится в Vportion, 
 и столько же эритроцитов попадает в мембрану, "не уместившийся обьем" переходил на следующий шаг, 
 (т.е. Vportion на следующем шаге равнялся не нулю , а этому оставшемуся обьему)
 оставшийся обьем Vportion не влиет на остальные параметры.
 """ 
 h=0 # шаг нефильтруемых
 warning=0
 while i <= (Tsim - 1):
  
  if Vwater < Vwaterall: # выполняется пока обьем прошедшей воды меньше обьема всей воды
      VEXIT[i] = VEXIT[i-1] +  a*d**2/4*math.pi*dt*(Np - ae)  #считается  обьем всей смеси-жидкости на выходе мембраны (на каждом шаге)
      Vwater =   Vwater     +  a*d**2/4*math.pi*dt*(Np - ae)  # считается обьем воды прошедшей через мембрану
      PotokBufer[i]=+  (a*d**2/4*math.pi*(Np - ae))/dVwater
  else:          
    VEXIT[i] = VEXIT[i-1] # если вся вода прошла то значение не изменяется
    Vwater = Vwater # если вся вода прошла то значение не изменяется
  if Vwater>Vwaterall:
      Vwater=Vwaterall
      VEXIT[i] = VEXIT[i-1]
  if p <= round(Ne):# если кол-во прошедших эритроцитов меньше всех эритроцитов
   
  
 
    """
     здесь начинается стохастическое опредление закрытых пор
    """
    
  filled = St[St != 0] #занятые поры
  empty = Np - filled.shape[0] #считаем, сколько пор свободно
  if Vwater<Vwaterall:
   Esred=round(Gemotocrid/100*(VEXIT[i]-VEXIT[i-1])/Ve)
   if AEM[i-1]==0:
     Esred=round(Gemotocrid/100*(VEXIT[i]-VEXIT[i-1]+Voverpor)/Ve)

   if empty ==0:
        Varaity=0 
   else:
     Varaity=(Esred/empty)
    
   close=black_bars2(Varaity, empty) #Теперь это не список закрытых а их число 
  else: 
         close=npe
         if npe>empty:
             close=empty
  Newempty=np.zeros((empty-close))
  AE=close
  if AE > (Ne - ae - p):
      AE = Ne - p - ae
  if AE==0:
    Voverpor=Voverpor+VEXIT[i]-VEXIT[i-1]
  if AE!=0:
    Voverpor=0
  if empty ==0:
      Varaity=0 
  NeededAmount=close
  if npe != 0:
   Teapr, N1gone, N2gone = Monte_Curve(te, te1, te2, npe, N1, N2, percentbadf, NeededAmount)
   if N2gone>N2:
       N2gone=N2
   if N1gone>N1:
        N1gone=N1   
   BadErinMembr[i]=N1gone
   GoodErinMembr[i]=N2gone
   AE = len(Teapr)
  else:
      Teapr=np.zeros(NeededAmount)
  Newclose=np.concatenate((Teapr, Newempty)) 
    
  npe=npe-AE
  N1=N1-N1gone
  N2=N2-N2gone
    
  St = np.concatenate((filled, Newclose)) 

  ae=ae+AE
  if ae>Np:
      ae=Np
    

  St[St> 0] = St[St > 0]-dt # изначально не было прибавки времени для закрытх пор


  passed = St[St < 0].shape[0]
  if passed > ae:
      passed=ae
  p += passed
  ae -= passed
  VEXIT[i] += passed*Ve
  PotokEr[i] += passed*Ve/dVwater/dt
  St[St < 0] = 0
  dV[i] = (VEXIT[i] - VEXIT[i - 1])/dt/dVwater #Скорсть течения жидкости
  AEM[i] = ae   # считается количество эритроцитов в мембране (на каждом шаге)
  Gemotocrit[i]=Gemotocrid
  AET[i]=AE
  t[i] = t[i - 1] + dt # считаетя время на следующем шаге
  if npe !=0:
   n2=N2/npe*100
  else:
     n2=0 
  """вывод основных параметров, для отслеживания работы программы"""
  if i % 1==0:  
   X = str( "Начальный Гематокрит"+str(Gemostart)+"%"+'\n'\
         +"Время " +  str(t[i]) + '\n' + "Отношение вышедшего обьема к входному в % "\
        + str(round(VEXIT[i]/V*100, 4)) + '\n'+"Скорость течения жидкости "+str(round(dV[i], 5))+'\n'\
        +"Количество клеток в мембране в данное время через  ae"+ str(ae)+'\n'\
            +"Количество клеток в мембране в данное время через массив St "+ str(St[St != 0].shape[0])+'\n'\
            +"Количество новых клеток в мембране"+str(AE) +'\n'\
            +"Количество новых клеток 1го типа"+str(N1gone) +'\n'\
            +"Количество новых клеток 2го типа"+str(N2gone) +'\n'\
            +"Отношение плохофильтр. к пропускаемым"+str(N2gone/(N2gone+N1gone)*100) +'\n'\
            
            +"Количество ожидаемых клеток в мембране"+str(Esred) +'\n'\
            +"Количество новых пропущенных всего"+str(passed) +'\n'\
           
            +"Количество всего пропущенных клеток "+str(p) +  '\n'\
            +"Гематокрит прошедшей смеси в %"+str(passed*Ve/VEXIT[i]*100) +  '\n'\
            +"Отношения плохо фильтруемых к общему числу оставшихся в %" +str(n2) +  '\n'\
        + "i = " + str(i) + '\n'\
        + "Проверка постоянства обьема " + str(round((Ve*(Ne - p + ae) + VEXIT[i] + (Vwaterall - Vwater))/V, 5))+ '\n'\
            + "Гематокрит (через VEXIT[i]) "+str(round(GemoVEXIT, 5))) 
   label.config(text = str(X))
   window.update()
  """ тут вывод заканчивается"""
  #пересчитывем новый гематокрит 
  if  (V-Vwater-Ve*(npe)) >0:
       
   Gemotocrid = (npe)*Ve/(V-Vwater-Ve*(npe))*100
  else: Gemotocrid= 100
  #if Gemotocrid < 0: #Данные условия не требуются, из-за того что гематокрит считается другим способом, более точным
      #Gemotocrid=0
  #if Gemotocrid > 100: 
      #Gemotocrid=100 
  #Gemotocrid= GemoVEXIT
  if VEXIT[i]/V >= 0.99:#проверка не превысил ли вышедший обьем входной
   if ae==0:
       t[i] = t[i - 1] + dt
       VEXIT[i] = V
       print("Поступившая жидкость протекла через мембрану")
       break;
  if npe<=0:
      if ae==0:
          print("Поступившая жидкость протекла через мембрану")
          break;
       
  i = i + 1
 print("Симуляция закончилась на итерации №", i)
 print("Время за которое прошла смесь в сек", t[i])
 print("Проверка постоянства обьема " , str((Ve*(Ne - p + ae) + VEXIT[i] + (Vwaterall - Vwater))/V))  
 if Ne != 0:
     print("эритроцитов пропущено", p, " из ожидавщихся ", Ne, '\n'\
       ,"Процент прошедших ", p/Ne*100, "%") 
 print("Гематокрит на выходе", str(p*Ve/(VEXIT[i])*100))
 print("   ")
 Fname=name+".txt"
 F = open(Fname, 'w')
 elem="Гематокрит начальный "+str(Gemostart)+'\n'+"Симуляция закончилась на итерации №"+str(i)+'\n'\
     +"Время за которое прошла смесь в сек "+str(t[i])+'\n'\
     +"Гематокрит на выходе "+   str(p*Ve/(VEXIT[i])*100) +'\n'
 F.write(elem) 
 F.close()
 """далее идет блок, который удаляет не нужные элементы массивов, которые не несут информации"""
 VEXIT1 = []
 t1     = []
 dV1    = []

 VEXIT2 = []
 t2     = []
 dV2    = []
 Gemotocrit2=[]
 PotokBufer2=[]
 PotokEr2=[]
 AET2=[]
 AEM2=[]
 GoodErinMembr2=[]
 BadErinMembr2=[]

 AEM1   = []
 Gemotocrit1=[]
 AET1=[]
 dVc1=[]
 Gemotocritc1=[]
 tdV1=[]
 tGemo1=[]
 ij=0
 tmax=max(t)
 tmin=min(t)
 dVmax=max(dV)
 dVmin=min(dV[dV != 0])
 dVmin11=1.1*min(dV[dV != 0])
 Gmax=max(Gemotocrit)
 Gmin=min(Gemotocrit[Gemotocrit != 0])
 Gmin11=1.01 *min(Gemotocrit[Gemotocrit != 0])
 n=1

 VEXIT2.append(VEXIT[1])
 t2.append(t[1])
 dV2.append(dV[1])
 Gemotocrit2.append(Gemotocrit[1])
 PotokBufer2.append(PotokBufer[1])
 PotokEr2.append(PotokEr[1])
 AET2.append(AET[1])
 AEM2.append(AEM[1])
 GoodErinMembr2.append(GoodErinMembr[1])
 BadErinMembr2.append(BadErinMembr[1])

 VEXIT2.append(VEXIT[2])
 t2.append(t[2])
 dV2.append(dV[2])
 Gemotocrit2.append(Gemotocrit[2])
 PotokBufer2.append(PotokBufer[2])
 PotokEr2.append(PotokEr[2])
 AET2.append(AET[2])
 AEM2.append(AEM[2])
 GoodErinMembr2.append(GoodErinMembr[2])
 BadErinMembr2.append(BadErinMembr[2])


 VEXIT2.append(VEXIT[3])
 t2.append(t[3])
 dV2.append(dV[3])
 Gemotocrit2.append(Gemotocrit[3])
 PotokBufer2.append(PotokBufer[3])
 PotokEr2.append(PotokEr[3])
 AET2.append(AET[3])
 AEM2.append(AEM[3])
 GoodErinMembr2.append(GoodErinMembr[3])
 BadErinMembr2.append(BadErinMembr[3])

 VEXIT2.append(VEXIT[4])
 t2.append(t[4])
 dV2.append(dV[4])
 Gemotocrit2.append(Gemotocrit[4])
 PotokBufer2.append(PotokBufer[4])
 PotokEr2.append(PotokEr[4])
 AET2.append(AET[4])
 AEM2.append(AEM[4])
 GoodErinMembr2.append(GoodErinMembr[4])
 BadErinMembr2.append(BadErinMembr[4])

 VEXIT2.append(VEXIT[5])
 t2.append(t[5])
 dV2.append(dV[5])
 Gemotocrit2.append(Gemotocrit[5])
 PotokBufer2.append(PotokBufer[5])
 PotokEr2.append(PotokEr[5])
 AET2.append(AET[5])
 AEM2.append(AEM[5])
 GoodErinMembr2.append(GoodErinMembr[5])
 BadErinMembr2.append(BadErinMembr[5])

 VEXIT2.append(VEXIT[6])
 t2.append(t[6])
 dV2.append(dV[6])
 Gemotocrit2.append(Gemotocrit[6])
 PotokBufer2.append(PotokBufer[6])
 PotokEr2.append(PotokEr[6])
 AET2.append(AET[6])
 AEM2.append(AEM[6])
 GoodErinMembr2.append(GoodErinMembr[6])
 BadErinMembr2.append(BadErinMembr[6])

 VEXIT2.append(VEXIT[8])
 t2.append(t[8])
 dV2.append(dV[8])
 Gemotocrit2.append(Gemotocrit[8])
 PotokBufer2.append(PotokBufer[8])
 PotokEr2.append(PotokEr[8])
 AET2.append(AET[8])
 AEM2.append(AEM[8])
 GoodErinMembr2.append(GoodErinMembr[8])
 BadErinMembr2.append(BadErinMembr[8])

 while ij<=i-1:  
  VEXIT1.append(VEXIT[ij])
  t1.append(t[ij])
  dV1.append(dV[ij])
  AEM1.append(AEM[ij])
  Gemotocrit1.append(Gemotocrit[ij])
  AET1.append(AET[ij])
  if Gemotocrit[ij] != 0:
     if Gemotocrit[ij] != Gmax:
      if Gemotocrit[ij] != Gmin:
        if Gemotocrit[ij] > Gmin11:
         tGemo1.append(t[ij])
         Gemotocritc1.append(Gemotocrit[ij]) 
  if dV[ij] != 0:
     if dV[ij] != dVmax:
      if dV[ij] != dVmin:
        if dV[ij] > dVmin11:
         if t[ij]>tmin*0.02:
          if t[ij]<= tmax*0.98:
           tdV1.append(t[ij])
           dVc1.append(dV[ij])
  if ij // 5 ==n:
     n+=1
     VEXIT2.append(VEXIT[ij])
     t2.append(t[ij])
     dV2.append(dV[ij])
     Gemotocrit2.append(Gemotocrit[ij])
     PotokBufer2.append(PotokBufer[ij])
     PotokEr2.append(PotokEr[ij])
     AET2.append(AET[ij])
     AEM2.append(AEM[ij])
     GoodErinMembr2.append(GoodErinMembr[ij])
     BadErinMembr2.append(BadErinMembr[ij])
  
  ij +=1

 return VEXIT1, t1, dV1, AEM1, Gemotocrit1, AET1, tdV1, dVc1, tGemo1, Gemotocritc1, VEXIT2, t2, dV2, Gemotocrit2, PotokBufer2, PotokEr2, AET2, AEM2, GoodErinMembr2, BadErinMembr2

def SimForStSpeed(V, Np, Gemotocrid, percent, percentbadf, Percentdouble,  Tsim, dt, a , te, te2, name, tcontrol):
 dVwater =a*d**2/4*math.pi*dt*(Np)/dt
 Ve = 100*10**(-18) *10**9#обьем эритроцита
 Ne = int(round(Gemotocrid/100*V/Ve))
 print("Гематокрит", Gemotocrid)
 print("время прохождения эритроцитом поры", te)
 print("время прохождения плохофильтр. эритроцитом поры", te2)
 Percentdouble=Percentdouble/100
 DoubleNp=round(Percentdouble*Np)

 
 # так тоже читаемо
 VEXIT, dV, t, AEM, AET ,Gemotocrit, PotokBufer, PotokEr= (np.zeros(Tsim) for _ in range(8))
 nfe=int(percent/100*Ne) # число нефильтруемых
 print("отношение нефильтруемых к порам", nfe/Np)
 nfeNum = np.linspace(1, int(Ne), int(nfe))
 for i in range(len((nfeNum))):
    nfeNum[i]=round(nfeNum[i])
 pae=0 # общее число попавших в мембрану
 npe=Ne # еще не прошедшие через мембрану клетки
 p = 0 # количество уже прошедших эритроцитов через мембрану
 ae = 0 # количество активных эритроцитов в мембране 
 blockae=0
 blockperc=0
 St  = np.zeros(Np) #массив отвеающий за состояние каждой поры (открыта - 0, занята - любое отличное от нуля число) 
 i = 1 # начальный шаг
 Voverpor=0 
 Vwater = 0 # обьем воды прошедшей через мембрану
 Vwaterall = V - Ne*Ve # обьем всей воды участвующей в расчетах, в не зависимости от состояния
 GemoVEXIT=Gemotocrid
 Gemostart=Gemotocrid
 N1=Ne*(1-percentbadf/100)
 N2=Ne*(percentbadf/100)
 """
 я высчитывал пройдет или непройдет эритроцит через мембрану при помощи таких этих порций
 суть в следующем, если данный обьем порции кратен целому числу обьемов воды,
 содержащих только 1 эритроцит (это число я считал ниже исходя из равномерной распределенности эритроцитов в жидкости)
 то считается сколько таких целых обьемов уместится в Vportion, 
 и столько же эритроцитов попадает в мембрану, "не уместившийся обьем" переходил на следующий шаг, 
 (т.е. Vportion на следующем шаге равнялся не нулю , а этому оставшемуся обьему)
 оставшийся обьем Vportion не влиет на остальные параметры.
 """ 
 h=0 # шаг нефильтруемых
 warning=0
 while i <= (Tsim - 1):
  
  if Vwater < Vwaterall: # выполняется пока обьем прошедшей воды меньше обьема всей воды
      VEXIT[i] = VEXIT[i-1] +  a*d**2/4*math.pi*dt*(Np - ae)  #считается  обьем всей смеси-жидкости на выходе мембраны (на каждом шаге)
      Vwater =   Vwater     +  a*d**2/4*math.pi*dt*(Np - ae)  # считается обьем воды прошедшей через мембрану
      PotokBufer[i]=+  (a*d**2/4*math.pi*(Np - ae))/dVwater
  else:          
    VEXIT[i] = VEXIT[i-1] # если вся вода прошла то значение не изменяется
    Vwater = Vwater # если вся вода прошла то значение не изменяется
  if Vwater>Vwaterall:
      Vwater=Vwaterall
      VEXIT[i] = VEXIT[i-1]
  if p <= round(Ne):# если кол-во прошедших эритроцитов меньше всех эритроцитов
   
  
 
    """
     здесь начинается стохастическое опредление закрытых пор
    """
         
    filled = St[St != 0] #занятые поры
    empty = Np - filled.shape[0] #считаем, сколько пор свободно

    Esred=round(Gemotocrid/100*(VEXIT[i]-VEXIT[i-1])/Ve)
    if AEM[i-1]==0:
       Esred=round(Gemotocrid/100*(VEXIT[i]-VEXIT[i-1]+Voverpor)/Ve)
    if Vwater>=Vwater:
           AE=npe
           if npe>empty:
               AE=empty
    if empty ==0:
          Varaity=0 
    else:
       Varaity=(Esred/empty)
      
    close=black_bars2(Varaity, empty) #Теперь это не список закрытых а их число 
    Newempty=np.zeros((empty-close))
    AE = close
    NeededAmount=AE
    Teapr, dt2, N1gone, N2gone = Monte_Curve(te, te1, te2, npe, N1, N2, percentbadf, NeededAmount, dt)
      
    Newclose=np.concatenate((Teapr, Newempty)) 
      
    npe=npe-AE
    N1=N1-N1gone
    N2=N2-N2gone
      
    St = np.concatenate((filled, Newclose)) 

    ae=ae+AE
      

    St[St> 0] = St[St > 0]-dt # изначально не было прибавки времени для закрытх пор


    passed = St[St < 0].shape[0]
    if passed > ae:
        passed=ae
    p += passed
    ae -= passed
    VEXIT[i] += passed*Ve
    PotokEr[i] += passed*Ve/dVwater/dt
    St[St < 0] = 0
    dV[i] = (VEXIT[i] - VEXIT[i - 1])/dt/dVwater #Скорсть течения жидкости
    AEM[i] = ae   # считается количество эритроцитов в мембране (на каждом шаге)
    Gemotocrit[i]=Gemotocrid
    AET[i]=AE
    t[i] = t[i - 1] + dt # считаетя время на следующем шаге
    """вывод основных параметров, для отслеживания работы программы"""
    if i % 1==0:  
      X = str( "Начальный Гематокрит"+str(Gemostart)+"%"+'\n'\
            +"Время " +  str(t[i]) + '\n' + "Отношение вышедшего обьема к входному в % "\
           + str(round(VEXIT[i]/V*100, 4)) + '\n'+"Скорость течения жидкости "+str(round(dV[i], 5))+'\n'\
           +"Количество клеток в мембране в данное время через  ae"+ str(ae)+'\n'\
               +"Количество клеток в мембране в данное время через массив St "+ str(St[St != 0].shape[0])+'\n'\
               +"Количество новых клеток в мембране"+str(AE) +'\n'\
               +"Количество новых клеток 1го типа"+str(N1gone) +'\n'\
               +"Количество новых клеток 2го типа"+str(N2gone) +'\n'\
                
               
               +"Количество ожидаемых клеток в мембране"+str(Esred) +'\n'\
               +"Количество новых пропущенных всего"+str(passed) +'\n'\
              
               +"Количество всего пропущенных клеток "+str(p) +  '\n'\
               +"Гематокрит прошедшей смеси в %"+str(passed*Ve/VEXIT[i]*100) +  '\n'\

           + "i = " + str(i) + '\n'\
           + "Проверка постоянства обьема " + str(round((Ve*(Ne - p + ae) + VEXIT[i] + (Vwaterall - Vwater))/V, 5))+ '\n'\
               + "Гематокрит (через VEXIT[i]) "+str(round(GemoVEXIT, 5))) 
    label.config(text = str(X))
    window.update()
    """ тут вывод заканчивается"""
    #пересчитывем новый гематокрит 

    GemoVEXIT=(npe)*Ve/(V-VEXIT[i]-ae*Ve)*100
    Gemotocrid=GemoVEXIT
    Gemotocrid = (npe)*Ve/(V-Vwater-Ve*(p+ae))*100
  if t[i]>=tcontrol:#проверка не превысил ли вышедший обьем входной
   Vstationar=dV[i]
   print("Стационарная скорость установилась в момент времени", tcontrol)
   break;
  i=i+1       
 return Vstationar

"""функция котроая сохраняет все данные в виде txt"""
def save(t2, VEXIT2, dV2, Gemotocrit2, PotokBufer2, PotokEr2, AET2, AEM2, GoodErinMembr2, BadErinMembr2, name):      
   name2=str(name)+'.xlsx'
   df = ({'Время в сек.': list(t2),
                   'Объем в мкл.': list(VEXIT2),
                   'Относительная скорость течения': list(dV2),
                   'Гематокрит в %': list(Gemotocrit2),
                   'Поток буфера в отн. ед.': list(PotokBufer2),
                   'Поток клеток в отн. ед.': list(PotokEr2),
                   'Количесвто поступивших клеток в мембрану': list(AET2),
                   'Количество  хороших клеток в мембране': list(GoodErinMembr2),
                   'Количество плохих клеток в мембране': list(BadErinMembr2)})
                   
   df = pd.DataFrame(data=df)
   writer=pd.ExcelWriter(name2, engine='xlsxwriter')
   df.to_excel(writer, sheet_name='Данные', index=False)
   writer.close()

def savestatspeed(tefiltr, dVstat, name):      
   name2=str(name)+'.xlsx'
   df = ({'Время фильтрации эритроцита в сек.': list(tefiltr),
                   'Стационарная скорость': list(dVstat)})
                   
   df = pd.DataFrame(data=df)
   writer=pd.ExcelWriter(name2, engine='xlsxwriter')
   df.to_excel(writer, sheet_name='Данные', index=False)
   writer.close()

 
def graphs(VEXIT1, t1, dV1, AEM1, Gemotocrit1, AET1, tdV1, dVc1, tGemo1, Gemotocritc1, name):

 """построение основных графиков"""
 plt.plot(t1, dV1,'.')
 plt.legend(title=name)
 plt.grid()
 plt.xlabel("Время в с",fontsize = 10)
 plt.ylabel("скорость течения жидкости в м^(3)/сек", fontsize = 10 )
 plt.title("скорость течения жидкости", fontsize = 10)
 NAME=str(name)+' скорость течения жидкости'+'.png'
 plt.savefig(NAME, dpi = 200)
 plt.show()
 
 
 plt.plot(t1, VEXIT1,'.') 
 plt.legend(title=name)
 plt.grid()
 plt.xlabel("Время в с",fontsize = 10)
 plt.ylabel("Обьем жидкости регистрируемый за мембраной в м^3")
 plt.title("Обьем крови на выходе мембраны от времени", fontsize = 10)
 NAME=str(name)+' Обьем на выходе'+'.png'
 plt.savefig(NAME, dpi = 200)
 plt.show()

 
 plt.plot(t1, AEM1,'.') 
 plt.legend(title=name)
 plt.grid()
 plt.xlabel("Время в с",fontsize = 10)
 plt.ylabel("Количество Эритроцитов в мембране")
 plt.title("Количество Эритроцитов в мембране от времени", fontsize = 10)
 NAME=str(name)+' Количество Эритроцитов в мембране'+'.png'
 plt.savefig(NAME, dpi = 200)
 plt.show()
 
 plt.plot(t1, Gemotocrit1,'.') 
 plt.legend(title=name)
 plt.grid()
 plt.xlabel("Время в с",fontsize = 10)
 plt.ylabel("Гематокрит в %", fontsize = 10)
 plt.title("Гематокрит от времени", fontsize = 10)
 NAME=str(name)+' Гематокрит в процентах'+'.png'
 plt.savefig(NAME, dpi = 200)
 plt.show()
 
 plt.plot(t1, AET1,'.') 
 plt.legend(title=name)
 plt.grid()
 plt.xlabel("Время в с",fontsize = 10)
 plt.ylabel("Количество новых поступивших клеток в мембрану ", fontsize = 10)
 plt.title("Количество новых поступивших клеток в мембрану", fontsize = 10)
 NAME=str(name)+' Количество новых поступивших клеток в мембрану'+'.png'
 plt.savefig(NAME, dpi = 200)
 plt.show() 

 plt.plot(tdV1, dVc1,'.', color="green") 
 plt.legend(title=name)
 plt.grid()
 plt.xlabel("Время в с",fontsize = 10)
 plt.ylabel("скорость течения жидкости в м^(3)/сек", fontsize = 10 )
 plt.title("скорость течения жидкости", fontsize = 10)
 NAME=str(name)+' скорость течения жидкости (без 0)'+'.png'
 plt.savefig(NAME, dpi = 200)
 plt.show()
 
 plt.plot(tGemo1, Gemotocritc1,'.', color="green") 
 plt.legend(title=name)
 plt.grid()
 plt.xlabel("Время в с",fontsize = 10)
 plt.ylabel("Гематокрит", fontsize = 10)
 plt.title("Гематокрит от времени", fontsize = 10)
 NAME=str(name)+' Гематокрит (без 0)'+'.png'
 plt.savefig(NAME, dpi = 200)
 plt.show()
 
 """
 #разбиение графиков при помощи функции timesplit
 timesplit(VEXIT1, t1, "V регистрируемый секунд за мембраной")
 timesplit(dV1c, tdV1, "скорость течения жидкости")
 timesplit(AEM1, t1, "Количество Эритроцитов в мембране")
 timesplit(Gemotocrit1c, tGemo1, "Гемотокрит")   
 """
"""входные параметры"""
R=10**1

Tsim = 10**7
V = 2500*10**(-9)*10**9/R
Ve = 100*10**(-18)*10**9

Np = int(int(10**5)*6/R)
L = 4*10**(-6)*10**3 #длина поры
d = 3.5*10**(-6)*10**3 #диаметр поры
a = 1/1.41374728*10**3/32 # относительный параметр, отвечающий за скорость течения жидкости через одну пору
b = a/100*8.55 #относительный параметр, отвечающий за скорость продвижения эритроцита через мембрану (относительно скорости воды "а")

percentnonf=0#Процент нефильтруемых


Percentdouble=0
Gemotocrid=1

te = (L + 8*Ve/(math.pi*d**2))/b
 # время за которое эритроцит проходит через пору
te=round(te, 6)
te1=te
te2=te*100
dt = 10**(-1)*te



percentbadf=1
GetSImDATA(V, Np, Gemotocrid, percentnonf, percentbadf, Percentdouble,  Tsim, dt, a ,te, te2)

