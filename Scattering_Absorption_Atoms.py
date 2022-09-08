import random
import numpy as np
import matplotlib.pyplot as plot
import math
import sympy as sp

def stf(f,a,b):
    return eval(f,{'x':a,'y':b,'log10':np.log10})

pi=math.pi
n=3.5*10**(2)
N=10**6
k=8.617*10**(-5)
T=3600
r_H=120*10**(-10)
a0=5.29*10**(-9)

E_e_eV=5.11*10**5   
E_pr_eV=9.38*10**(8)
e=4.8*10.**(-10.)
c=3.*10.**(10.)
m_e=9.1*10.**(-28.)
r0=e**2/(m_e*c**2)
sigma_T=8./3.*pi*r0**2
r_proton=0.84*10**(-13)
sigma_proton=pi*r_proton**2
sigma=pi*a0**2

print("Set the population's spacial limits special with f(X,Y)=const1 as the upper limit and g(X,Y)=const2 as the lower limit.\n")

(x,y)=sp.symbols('x,y')

f=input("Import the equation f(x,y)= ")
Const1=float(input('Upper Limit Const1: '))
    
g=input("Import the equation g(x,y)= ")
Const2=float(input('Lower Limit Const2: '))    

x_spac=0*np.ones(shape=(N,1))
#y_spac=Const2*np.ones(shape=(N,1))
#y_spac=0.*np.ones(shape=(N,1))
y_spac=-Const1*np.ones(shape=(N,1))
theta=pi/2.*np.ones(shape=(N,1))

E0=1.54

E_ph=E0*np.ones(shape=(N,1))

Prob_H=[]
E_el=[]
E_pr=[]
H_energ=[]
min_tr=0.
for i in range(1,11):
    H_energ.append(13.6/i**2)
    Prob_H.append(np.exp(-(H_energ[0]-H_energ[i-1])/(k*T)))
    if 13.6/i**2-E0<=0.:
        E_el.append( E_e_eV+(E0-13.6/i**2)/2.)
        E_pr.append( E_pr_eV+(E0-13.6/i**2)/2.)
        if min_tr==0.:
            min_tr=i

Prob_H=Prob_H/sum(Prob_H)

n_atom=n*(sum(Prob_H[0:min_tr]))
n_electron=n*(1-sum(Prob_H[0:min_tr]))
n_proton=n_electron

l_atom=(n_atom*sigma)**(-1)/10**(11)
l_electron=(n_electron*sigma_T)**(-1)/10**(11)
l_proton=(n_proton*sigma_proton)**(-1)/10**(11)

col_number=np.zeros(shape=(N))

l_tot=(1./l_atom+1./l_electron+1./l_proton)**(-1.)
   
for i in range(0,N):
    if x_spac[i]!=None and y_spac[i]!=None:
        while stf(f,x_spac[i],y_spac[i])<=Const1 and stf(g,x_spac[i],y_spac[i])>=Const2:
            r=(l_atom**(-1)+l_electron**(-1)+l_proton**(-1))*random.random()
            if r<=(l_atom)**(-1):
                l=l_atom
            if r>(l_atom)**(-1) and r<=(l_electron)**(-1)+(l_atom)**(-1):
                l=l_electron
            if r>(l_electron)**(-1)+(l_atom)**(-1):
                l=l_proton
            r_collision=-l*(math.log(random.random()))
            y_spac[i]=y_spac[i]+r_collision*(math.sin(theta[i]))
            x_spac[i]=x_spac[i]+r_collision*(math.cos(theta[i]))
            if stf(f,x_spac[i],y_spac[i])<=Const1 and stf(g,x_spac[i],y_spac[i])>=Const2:
                if l==l_electron:                    
                    r=random.randint(0, len(E_el))
                    E_e_temp=E_el[r]
                    gamma_e=E_e_temp/E_e_eV
                    beta_e=np.sqrt(1-1/gamma_e**2)                    
                    theta_temp=2*math.pi*random.random()
                    theta_sc=np.abs(theta_temp-theta[i])                
                    E_ph[i]=E_ph[i]*(1-beta_e*np.cos(theta[i]))/(1-beta_e*np.cos(theta_temp)+E_ph[i]/E_e_temp*(1-np.cos(theta_sc)))
                    theta[i]=theta_temp
                if l==l_proton:
                    r=random.randint(0, len(E_el))
                    E_pr_temp=E_pr[r]
                    gamma_e=E_e_temp/E_pr_eV
                    beta_e=np.sqrt(1-1/gamma_e**2)                    
                    theta_temp=2*math.pi*random.random()
                    theta_sc=np.abs(theta_temp-theta[i])                
                    E_ph[i]=E_ph[i]*(1-beta_e*np.cos(theta[i]))/(1-beta_e*np.cos(theta_temp)+E_ph[i]/E_pr_temp*(1-np.cos(theta_sc)))
                    theta[i]=theta_temp
                else:
                    theta[i]=2*math.pi*random.random()                
                col_number[i]=col_number[i]+1

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
print("Theoretical percentage that don't collide: "+str(np.exp(-(Const1-Const2)/l_tot))) 
print("Accelerated particles percentage: "+str(np.count_nonzero(E_ph>E0)/N)) 
print("Absorbed particles percentage: "+str(np.count_nonzero(np.isnan(theta))/N)) 
print("Accelerated particles that got through percentage: "+str(len(half_acc)/N)) 
      
Angl=plot.figure()
plot.hist(theta,100, edgecolor='black', linewidth=1.2,color='white')
plot.xlabel('Theta(rad)')
plot.ylabel('N')
plot.savefig('PowerLaw_Theta_τ_'+str((Const1-Const2)/l_tot)+'.pdf')
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
plot.savefig('PowerLaw_Energ_τ_'+str((Const1-Const2)/l_tot)+'.pdf')
Energ=plot.figure()
plot.yscale('log')
plot.hist(half_energ,bins, edgecolor='black', linewidth=1.2,color='blue',alpha=0.8)
plot.hist(half_acc,bins, edgecolor='black', linewidth=1.2,color='magenta',alpha=0.7)
plot.hist(meas_energ,bins, edgecolor='black', linewidth=1.2,color='red',alpha=1.)
plot.xlabel('Energy(eV)')
plot.ylabel('log(N)')
plot.tight_layout()
plot.legend(['Radiation in θ>π','Accelerated Photons','Radiation in Δθ='+str(2*dw)+'rad'])
plot.savefig('PowerLaw_Energ_through_τ_'+str((Const1-Const2)/l_tot)+'.pdf')
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
plot.savefig('PowerLaw_logEnerg_τ_'+str((Const1-Const2)/l_tot)+'.pdf')
Half_Energ=plot.figure()
plot.yscale('log')
plot.hist(np.log10(half_energ),bins, edgecolor='black', linewidth=1.2,color='blue',alpha=0.8)
plot.hist(np.log10(half_acc),bins, edgecolor='black', linewidth=1.2,color='magenta',alpha=0.7)
plot.hist(np.log10(meas_energ),bins, edgecolor='black', linewidth=1.2,color='red',alpha=1.)
plot.xlabel('log10[Energy(eV)]')
plot.ylabel('log(N)')
plot.legend(['Radiation in θ>π','Accelerated Photons','Radiation in Δθ='+str(2*dw)+'rad'])
plot.tight_layout()
plot.savefig('PowerLaw_Half_Energ_τ_'+str((Const1-Const2)/l_tot)+'.pdf')