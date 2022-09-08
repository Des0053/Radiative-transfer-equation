import random
import numpy as np
import matplotlib.pyplot as plot
import math
import sympy as sp

def stf(f,a,b):
    return eval(f,{'x':a,'y':b,'log10':np.log10})


N=10**6

pi=math.pi
n=10**6
E_e_eV=5.11*10**5   
e=4.8*10.**(-10.)
c=3.*10.**(10.)
m_e=9.1*10.**(-28.)
r0=e**2/(m_e*c**2)
sigma_T=8./3.*pi*r0**2
        
print("Set the population's spacial limits special with f(X,Y)=const1 as the upper limit and g(X,Y)=const2 as the lower limit.\n")

(x,y)=sp.symbols('x,y')

f=input("Import the equation f(x,y)= ")
Const1=float(input('Upper Limit Const1: '))
    
g=input("Import the equation g(x,y)= ")
Const2=float(input('Lower Limit Const2: '))    

x_spac=0*np.ones(shape=(N,1))
y_spac=np.sqrt(Const2**2)*np.ones(shape=(N,1))
#y_spac=-Const1*np.ones(shape=(N,1))
theta=pi/2.*np.ones(shape=(N,1))

Emin=10**2*E_e_eV
Emax=1.5*10**5*E_e_eV
N0=10**(10)
s=2.

E0=1.5*10**(5)

E_ph=E0*np.ones(shape=(N,1))


def N_e(E):
    if E>Emax or E<Emin:
        return 0
    else:
        return N0*(E/Emin)**(-s)

N_tot=(N_e(Emin)+N_e(Emax))*(Emax-Emin)/2.

def sigma(E):
    a=E/E_e_eV
    if a<10**(-3):
        return sigma_T
    else:
        return pi*r0**2/(2.*a**3)*(a**2+8.*a-a**2/(2.*a+1.)**2+2.*(a**2-2.*a-2.)*math.log(2*a+1.))
    #return pi*r0**2*8./3.

def l(E):
    #gamma=E_e/E_e_eV
    #if gamma>1.:
        #beta=np.sqrt(1-1/gamma**2)
        #return (gamma*n*sigma(E))**(-1)/(3*10**(18))*np.sqrt(1+(1./gamma**2-1.)*(beta+np.cos(Theta_sc))**2/(1.+beta*np.cos(Theta_sc))**2)
        #return (n*sigma(E))**(-1)/(3*10**(18))
    #else:
    return (n*sigma(E))**(-1)/(3*10**(18))
    
col_number=np.zeros(shape=(N))

abs_prob=0.

#for i in range(0,N):
   # r_collision=-l[E_ph[i]]*(math.log(random.random()))
   # y_spac[i]=y_spac[i]+r_collision*(math.sin(theta[i]))
   # x_spac[i]=x_spac[i]+r_collision*(math.cos(theta[i]))
   # theta_e=2*pi*random.random()
   # theta_sc=np.abs(theta_e-theta[i])
   # E_e=Emin+(Emax-Emin)*random.random()
   # if stf(f,x_spac[i],y_spac[i])<=Const1 and stf(g,x_spac[i],y_spac[i])>=Const2:
    #    col_number[i]=col_number[i]+1
    
for i in range(0,N):
    if x_spac[i]!=None and y_spac[i]!=None:
        while stf(f,x_spac[i],y_spac[i])<=Const1 and stf(g,x_spac[i],y_spac[i])>=Const2:
            r_collision=-l(E_ph[i])*(math.log(random.random()))
            y_spac[i]=y_spac[i]+r_collision*(math.sin(theta[i]))
            x_spac[i]=x_spac[i]+r_collision*(math.cos(theta[i]))
            if stf(f,x_spac[i],y_spac[i])<=Const1 and stf(g,x_spac[i],y_spac[i])>=Const2:
                prob=random.random()            
                if prob<=abs_prob:
                    theta[i]=None
                    y_spac[i]=None
                    x_spac[i]=None
                    col_number[i]=None
                    break
                else:             
                    col_number[i]=col_number[i]+1  
                    E_e_temp=(Emin**(-s+1.)-random.random()*(Emin**(-s+1.)-Emax**(-s+1.)))**(1./(1.-s))
                    gamma_e=E_e_temp/E_e_eV
                    beta_e=np.sqrt(1-1/gamma_e**2)                    
                    theta_temp=2*math.pi*random.random()
                    theta_sc=np.abs(theta_temp-theta[i])                
                    E_ph[i]=E_ph[i]*(1-beta_e*np.cos(theta[i]))/(1-beta_e*np.cos(theta_temp)+E_ph[i]/E_e_temp*(1-np.cos(theta_sc)))
                    theta[i]=theta_temp

    else:
        break
    print(str(i+1)+'/'+str(N))


half_beam=0
half_energ=[]
meas_beam=0
meas_energ=[]
meas_energ_tot=0.
dw=10**(-7)
half_acc=[]

for i in range(0,len(theta)):
    if theta[i]<=(1./2.+dw)*pi and theta[i]>=(1./2.-dw)*pi:
        meas_beam=meas_beam+1
        meas_energ.append(float(E_ph[i]))
        meas_energ_tot=meas_energ_tot+float(E_ph[i])
    if theta[i]<=pi:
        half_beam=half_beam+1
        half_energ.append(float(E_ph[i]))
        if E_ph[i]>E0:
            half_acc.append(float(E_ph[i]))

flux, Energy=np.histogram(E_ph,np.linspace(float(min(E_ph)), float(max(E_ph)),100))
#%%        
print('Photon percentage in Δθ='+str(2*dw)+'rad: '+str(meas_beam/N))  
print("Photon percentage that don't collide: "+str(np.count_nonzero(col_number==0)/N)) 
print('Photon percentage that got through: '+str(half_beam/N))  
print("Theoretical percentage that don't collide: "+str(np.exp(-(Const1-Const2)/l(E0)))) 
print("Accelerated particles percentage: "+str(np.count_nonzero(E_ph>E0)/N)) 
print("Absorbed particles percentage: "+str(np.count_nonzero(np.isnan(theta))/N)) 
print("Accelerated particles that got through percentage: "+str(len(half_acc)/N)) 
      
Angl=plot.figure()
plot.hist(theta,100, edgecolor='black', linewidth=1.2,color='white')
plot.xlabel('Theta(rad)')
plot.ylabel('N')
plot.savefig('PowerLaw_Theta_τ_'+str((Const1-Const2)/l(E0))+'.pdf')
Energ=plot.figure()
plot.yscale('log')
bins = np.linspace(float(min(E_ph)), float(max(E_ph)), 100)
plot.hist(E_ph,bins, edgecolor='black', linewidth=1.2,color='white')
plot.hist(half_energ,bins, edgecolor='black', linewidth=1.2,color='blue',alpha=0.8)
plot.hist(half_acc,bins, edgecolor='black', linewidth=1.2,color='magenta',alpha=0.7)
plot.hist(meas_energ,bins, edgecolor='black', linewidth=1.2,color='red',alpha=1.)
plot.plot(Energy[0:len(Energy)-1],flux,'k',label='_nolegend_')
plot.xlabel('Energy(eV)')
plot.ylabel('log(N)')
plot.tight_layout()
plot.legend(['Total radiation','Radiation in θ>π','Accelerated Photons','Radiation in Δθ='+str(2*dw)+'rad'])
plot.savefig('PowerLaw_Energ_τ_'+str((Const1-Const2)/l(E0))+'.pdf')
Energ=plot.figure()
plot.yscale('log')
plot.hist(half_energ,bins, edgecolor='black', linewidth=1.2,color='blue',alpha=0.8)
plot.hist(half_acc,bins, edgecolor='black', linewidth=1.2,color='magenta',alpha=0.7)
plot.hist(meas_energ,bins, edgecolor='black', linewidth=1.2,color='red',alpha=1.)
plot.xlabel('Energy(eV)')
plot.ylabel('log(N)')
plot.tight_layout()
plot.legend(['Radiation in θ>π','Accelerated Photons','Radiation in Δθ='+str(2*dw)+'rad'])
plot.savefig('PowerLaw_Energ_through_τ_'+str((Const1-Const2)/l(E0))+'.pdf')
Energ=plot.figure()
plot.yscale('log')
bins = np.linspace(float(min(np.log10(E_ph))), float(max(np.log10(E_ph))), 100)
plot.hist(np.log10(E_ph),bins, edgecolor='black', linewidth=1.2,color='white')
plot.hist(np.log10(half_energ),bins, edgecolor='black', linewidth=1.2,color='blue')
plot.hist(np.log10(half_acc),bins, edgecolor='black', linewidth=1.2,color='magenta',alpha=0.7)
plot.hist(np.log10(meas_energ),bins, edgecolor='black', linewidth=1.2,color='red')
plot.xlabel('log10[Energy(eV)]')
plot.ylabel('log(N)')
plot.tight_layout()
plot.legend(['Total radiation','Radiation in θ>π','Radiation in Δθ='+str(2*dw)+'rad','Accelerated Photons'])
plot.savefig('PowerLaw_logEnerg_τ_'+str((Const1-Const2)/l(E0))+'.pdf')
Half_Energ=plot.figure()
plot.yscale('log')
plot.hist(np.log10(half_energ),bins, edgecolor='black', linewidth=1.2,color='blue',alpha=0.8)
plot.hist(np.log10(half_acc),bins, edgecolor='black', linewidth=1.2,color='magenta',alpha=0.7)
plot.hist(np.log10(meas_energ),bins, edgecolor='black', linewidth=1.2,color='red',alpha=1.)
plot.xlabel('log10[Energy(eV)]')
plot.ylabel('log(N)')
plot.legend(['Radiation in θ>π','Accelerated Photons','Radiation in Δθ='+str(2*dw)+'rad'])
plot.tight_layout()
plot.savefig('PowerLaw_Half_Energ_τ_'+str((Const1-Const2)/l(E0))+'.pdf')