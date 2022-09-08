import random
import numpy as np
import matplotlib.pyplot as plot


Emin=10**6
Emax=10**(12)
N0=10**(20)
N=10**(6)
s=2.

def N_e(E):
    if E>Emax or E<Emin:
        return 0
    else:
        return N0*(E/Emin)**(-s)

#N_tot=(N_e(Emin)+N_e(Emax))*(Emax-Emin)/2.
Ep=np.linspace(Emin, Emax,10**3)
katan=[]
katan_numb=np.zeros(shape=(len(Ep),1))

for i in range(0,N):
    katan.append( (Emin**(-s+1.)-random.random()*(Emin**(-s+1.)-Emax**(-s+1.)))**(1./(1.-s)))
    for j in range(0,len(Ep)-1):
        if katan[i]>=Ep[j] and katan[i]<=Ep[j+1]:
            katan_numb[j]=katan_numb[j]+1
    print(str(i+1)+'/'+str(N))
#plot.hist(katan,bins=100)    
#plot.plot(Ep,N*(s-1.)/(Emin**(-s+1.)-Emax**(-s+1.))*Ep**(-s))

plot.plot(np.log10(Ep),np.log10(katan_numb))
plot.plot(np.log10(Ep),np.log10(N*(s-1.)/(Emin**(-s+1.)-Emax**(-s+1.))*Ep**(-s+1)),np.log10(Ep),np.log10(katan_numb))