import random
import numpy as np
import matplotlib.pyplot as plot
import math
import sympy as sp

def stf(f,a,b):
    return eval(f,{'x':a,'y':b,'log10':np.log10})

pi=math.pi
n=10**5
D=10.**3
N=10**5
k=8.617*10**(-5)
T=3600
r_H=120*10**(-12)
a0=5.29*10**(-9)

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
y_spac=-np.sqrt(Const1)*0*np.ones(shape=(N,1))
theta=pi/2.*np.ones(shape=(N,1))

Emin=10**0*E_e_eV
Emax=1.1*10**0*E_e_eV
N0=10**(10)

E0=10**7

E_ph=E0*np.ones(shape=(N,1))

def sigma(E):
    a=E/E_e_eV
    if a<10**(-3):
       return sigma_T
    else:
       return pi*r0**2/(2.*a**3)*(a**2+8.*a-a**2/(2.*a+1.)**2+2.*(a**2-2.*a-2.)*math.log(2*a+1.))
   
def l(E):
    return(n*sigma(E))**(-1)/(3*10**(18))

col_number=np.zeros(shape=(N))

abs_prob=0.0

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
                    if stf(f,x_spac[i],y_spac[i])<=Const1 and stf(g,x_spac[i],y_spac[i])>=Const2:
                        col_number[i]=col_number[i]+1
                        E_e_temp=-k*T*np.log(np.exp(-Emin/(k*T))-random.random()*(np.exp(-Emin/(k*T))-np.exp(-Emax/(k*T))))
                        gamma_e=E_e_temp/E_e_eV
                        beta_e=np.sqrt(1-1/gamma_e**2)
                        theta_temp=2*math.pi*random.random()
                        theta_sc=np.abs(theta_temp-theta[i])
                        E_ph[i]=E_ph[i]*(1-beta_e*np.cos(theta[i]))/(1-beta_e*np.cos(theta_temp)+E_ph[i]/E_e_temp*(1-np.cos(theta_sc)))
                        theta[i]=theta_temp

    else:
        break
    print(str(i+1)+'/'+str(N))

meas_beam=0
meas_energ=[]

for i in range(0,len(theta)):
    if theta[i]<=(1./2.+0.01)*pi and theta[i]>=(1./2.-0.01)*pi:
        meas_beam=meas_beam+1
        meas_energ.append(float(E_ph[i]))
        
print(meas_beam/N)  
      
Angl=plot.figure()
plot.hist(theta,100, edgecolor='black', linewidth=1.2,color='white')
plot.xlabel('Theta(rad)')
plot.savefig('MaxBol_Theta.pdf')
Energ=plot.figure()
plot.yscale('log')
plot.hist(E_ph,100, edgecolor='black', linewidth=1.2,color='white')
plot.hist(meas_energ,100, edgecolor='black', linewidth=1.2,color='blue')
plot.xlabel('Energy(eV)')
plot.savefig('MaxBol_Energ.pdf')
