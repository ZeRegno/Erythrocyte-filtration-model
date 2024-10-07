import numpy as np
import matplotlib.pyplot as plt
import math
from tkinter import Tk, Label
import pandas as pd




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
"""
Эта функция - основная. При помощи нее и ведется весь расчет
На вход подаются следующие значения:
    V - обьем всей жидкости - смеси
    Np - число пор
    Gemotocrid - значение гематокрита смеси
    Tsim - ограничивающие число "шагов" модели, по его достижение программа останавливается
      так же Tsim задает количество элментов всех массивов данных, 
      в последсвтии все не несущие информации элементы удаляются
    dt - шаг времени, от него зависит точность и время вычислений   
"""              
def Sim(V, Np, Gemotocrid, percent, percentbadf, Percentdouble,  Tsim, dt, a , te, te2, name):
 dVwater =a*d**2/4*math.pi*dt*(Np)/dt
 Ve = 100*10**(-18) *10**9#обьем эритроцита
 Ne = int(round(Gemotocrid/100*V/Ve))
 print("Гематокрит", Gemotocrid)
 print("время прохождения эритроцитом поры", te)
 print("время прохождения плохофильтр. эритроцитом поры", te2)
 Percentdouble=Percentdouble/100
 DoubleNp=round(Percentdouble*Np)

 
 # так тоже читаемо
 VEXIT, dV, t, AEM, AET ,Gemotocrit, PotokBufer, PotokEr, BadErinMembr, GoodErinMembr= (np.zeros(Tsim) for _ in range(10))
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
    
    # создание маасива двойных пор
    Doubleemptyamount=round(Percentdouble*empty)
    Doubleempty=np.zeros(empty)
    A=0
    while A<Doubleemptyamount:
        Doubleempty[A]=1
        A=A+1
    np.random.shuffle(Doubleempty)
    
    Esred1=round((Gemotocrid*(1-percentbadf/100))/100*(VEXIT[i]-VEXIT[i-1])/Ve)
    Esred=round(Gemotocrid/100*(VEXIT[i]-VEXIT[i-1])/Ve)
    if AEM[i-1]==0:
     Esred1=round((Gemotocrid*(1-percentbadf/100))/100*(VEXIT[i]-VEXIT[i-1]+Voverpor)/Ve)
    if Vwater>=Vwater:
        if VEXIT[i]-VEXIT[i-1] ==0:
         AE=npe
         if npe>empty:
             AE=empty
    if empty ==0:
        Varaity=0 
    else:
     Varaity=(Esred1/empty)
    
    #Пуассон, mu - средее число эритроцитов в мембране - вероятность закрытия поры
    #Вероятность закрытия должна меняться от концентрации, поэтому гематокрит надо пересчитывать!
    close=black_bars(Varaity, empty, dt)
    AE1 = close.nonzero()[0].shape[0] #число новых эритроцитов в мембране
    if percentbadf==0:
        Amountofbadf=0
        Varaity2=0
        
    else:
     Varaity2=(Esred/empty*percentbadf/100)
     Amountofbadf=black_bars2(Varaity2, empty) 
    AE=AE1+Amountofbadf
    BadErinMembr[i]=BadErinMembr[i-1]+Amountofbadf
    GoodErinMembr[i]=GoodErinMembr[i-1]+AE1
    g=0
    ij=0
    if Amountofbadf>0:
     for ij in range(len(close)):
      if close[ij]==0:
         close[ij]=-dt/2
         g=g+1
      if g>= Amountofbadf:
          break;
     np.random.shuffle(close)
        
      
    pae=pae+AE
    l=0
    Ll=0
    blockae=0
    blocked=St[St == -1].shape[0]
    if percent !=0:
     if Gemotocrid !=0:
      if blocked<(Np-DoubleNp):
       if h<=nfe:
        l=0
       if pae >= nfeNum[h]:
        while pae >= nfeNum[h]:
         h=h+1
         l=l+1
       X=0
       for ij in range(len(close)):
        if X<l:
         if close[ij]==dt/2 :
            close[ij]=-1
            if Doubleempty[ij]==1:
                close[ij]=dt/2
                Ll=Ll+1
            X=X+1 
     
     blockae=blockae+l-Ll
    
    npe=npe-AE
    #print(AE)
    St = np.concatenate((filled, close)) 
    #ae= St[St != 0].shape[0]
    ae=ae+AE
    if blockae > ae:
        print("ошибка заблокированных оказалось больше занятых!!!!!!")
    
  
  St[St>= dt] = St[St >=dt]+dt # изначально не было прибавки времени для закрытх пор
  St[St == dt/2] = dt # такая сложность нужна чтобы не прибавлять только что закрытым порам еще шаг времени
  St[St<= -dt]=St[St<= -dt]-dt
  St[St== -dt/2]=-dt
  
  passed1 = St[St >= te].shape[0]
  if passed1 > ae:
      passed1=ae
  p += passed1
  ae -= passed1
  VEXIT[i] += passed1*Ve
  PotokEr[i] += passed1*Ve/dVwater/dt
  St[St >= te] = 0
  
  passed2 = St[St <= -1*te2].shape[0]
  if passed2 > ae:
      passed2=ae
  p += passed2
  ae -= passed2
  VEXIT[i] += passed2*Ve
  PotokEr[i]+= passed2*Ve/dVwater/dt
  St[St <= -1*te2] = 0
  passed=passed1+passed2
  BadErinMembr[i]=BadErinMembr[i]-passed1
  GoodErinMembr[i]=GoodErinMembr[i]-passed2
  """
  passed=0
  passed2=0
  for ij in range(Np):
      if St[ij] != 0:
          if St[ij].real == dt:
              St[ij]=St[ij]+dt/2
          if St[ij].real == dt/2:
              St[ij]=St[ij]+dt/2
          if St[ij]>=te:
              if St[ij].imag==0:
                  St[ij]=0
                  passed=passed+1
          if St[ij]>=te2:
              if St[ij].imag==1:
                  St[ij]=0
                  passed=passed+1
                  passed2=passed2+1
  if passed>ae:
   passed=ae  
  p+=passed             
  ae-=passed
  VEXIT[i] += passed*Ve
  """
  if warning==0:
   if p>Ne:
      print("Ошибка, эритроцитов прошло больше чем есть")
      warning=1
  #print(VEXIT[i]/V0)
  blocked=St[St == -1].shape[0]
  
  if blocked>(Np-DoubleNp):
    ij=0
    while blocked>(Np-DoubleNp):
        if St[ij]==-1:
            St[ij]=0
            AE=AE+1
            passed=passed+1
            VEXIT[i] += Ve
            blocked=St[St == -1].shape[0]
        ij=ij+1
  dV[i] = (VEXIT[i] - VEXIT[i - 1])/dt/dVwater #Скорсть течения жидкости
  AEM[i] = ae   # считается количество эритроцитов в мембране (на каждом шаге)
  Gemotocrit[i]=Gemotocrid
  AET[i]=AE
  t[i] = t[i - 1] + dt # считаетя время на следующем шаге
  if p!=0:
      blockperc=round(blocked/(p+ae)*100,5)
  Newblper=0
  if AE!=0: 
    Newblper=round(l/AE*100, 3)
  Propush=0
  if Ne !=0:
      Propush=round(p/Ne*100,2)
  """вывод основных параметров, для отслеживания работы программы"""
  if i % 1==0:  
   X = str( "Начальный Гематокрит"+str(Gemostart)+"%"+'\n'\
         +"Время " +  str(t[i]) + '\n' + "Отношение вышедшего обьема к входному в % "\
        + str(round(VEXIT[i]/V*100, 4)) + '\n'+"Скорость течения жидкости "+str(round(dV[i], 5))+'\n'\
        +"Количество заблокированных пор "+str(blocked)+'\n'\
        +"Количество клеток в мембране в данное время через  ae"+ str(ae)+'\n'\
            +"Количество клеток в мембране в данное время через массив St "+ str(St[St != 0].shape[0])+'\n'\
            +"Количество новых клеток в мембране"+str(AE) +'\n'\
            +"Количество новых клеток 2го типа"+str(Amountofbadf) +'\n'\
             +"Количество новых заблокированных"+str(l-Ll) +'\n'\
            +"Количество Нефильтруемых попывшие в двойные"+str(Ll) +'\n'\
            +"Количество ожидаемых клеток в мембране"+str(Esred) +'\n'\
            +"Количество новых пропущенных всего"+str(passed) +'\n'\
            +"Количество новых пропущенных 1го типа"+str(passed1) +'\n'\
            +"Количество новых пропущенных 2го типа"+str(passed2) +'\n'\
            +"Вероятность попадания 2го типа"+str(Varaity2) +'\n'\
            +"Количество новых заблокированных в % от новых эритроцитов "+str(Newblper) +'\n'\
            +"Количество всего пропущенных клеток "+str(p) +  '\n'\
            +"Заблокированные в проценте от прошедших "+str(blockperc) +  '\n'\
            +"Количество всего пропущенных клеток в % "+str(Propush) +  '\n'\
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
  #if Gemotocrid < 0: #Данные условия не требуются, из-за того что гематокрит считается другим способом, более точным
      #Gemotocrid=0
  #if Gemotocrid > 100: 
      #Gemotocrid=100 
  #Gemotocrid= GemoVEXIT
  if VEXIT[i]/V >= 0.97:#проверка не превысил ли вышедший обьем входной
   if ae==0:
       t[i] = t[i - 1] + dt
       VEXIT[i] = V
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
Tsim = 10**7
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
Gemotocrid=1
percentbadf=0.1
Percentdouble=0

te = (L + 8*Ve/(math.pi*d**2))/b # время за которое эритроцит проходит через пору
te=round(te, 6)
te2=te*10
dt = 10**(-1)*te


percentbadf=0.1

GetSImDATA(V, Np, Gemotocrid, percentnonf, percentbadf, Percentdouble,  Tsim, dt, a ,te, te2)


