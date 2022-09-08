import random
import numpy as np
import matplotlib.pyplot as plot

c=3.*10.**(10.)
m_e=0.511*10**(6)/c**2
Emin=10**(-4)
Emax=10**(-1)
k=8.617*10**(-5)
T=300
N=10**(6)

Ep=np.linspace(Emin, Emax,10**3)
katan=[]
katan_numb=np.zeros(shape=(len(Ep),1))

for i in range(0,N):
    katan.append(-k*T*np.log(np.exp(-Emin/(k*T))-random.random()*(np.exp(-Emin/(k*T))-np.exp(-Emax/(k*T)))))
    for j in range(0,len(Ep)-1):
        if katan[i]>=Ep[j] and katan[i]<=Ep[j+1]:
            katan_numb[j]=katan_numb[j]+1
            break
    print(str(i+1)+'/'+str(N))

#plot.hist(katan,bins=100)    
#plot.plot(Ep,(s-1.)/(Emin**(-s+1.)-Emax**(-s+1.))*Ep**(-s))

plot.plot(np.log10(Ep),np.log10(katan_numb))