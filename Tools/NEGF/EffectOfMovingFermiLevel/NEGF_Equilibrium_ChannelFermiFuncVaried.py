#!/usr/bin/env python
# coding: utf-8

# In[286]:


import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
plt.rcParams.update({'font.size': 12})
import seaborn as sns
import cmath
import math 
def custom_plot_single(x_lst, y_lst, label_lst, xlim, ylim, label,pltname,
                       color=['k','r','b','g','o','br'],
                       linestyle=['solid','dashed','solid','dashed','solid','dashed'],
                       markertype=[None,None,'o','^','o','^'],
                       fillstyle=['none','none','none','none','full','full'],
                       linewidth=20*[3],
                       markevery=[45,30,50,40,56,72,63,95]):
    
    fig, ax1 = plt.subplots(1,1)
    fig.set_size_inches(10,10)
    fig.patch.set_facecolor('white')
    ax1.patch.set_facecolor('white')
    for p in range(0, len(x_lst)):
        ax1.plot(x_lst[p], y_lst[p], color[p],
                 linewidth=linewidth[p],
                 linestyle=linestyle[p],
                 marker=markertype[p],
                 fillstyle=fillstyle[p],
                 markevery=markevery[p],
                 markersize=8,
                 label=label_lst[p])

    ax1.legend(prop={'size': 20},loc='best')
    ax1.tick_params(which='minor', width=2, length=4, color='k')
    ax1.tick_params(which='major', width=2, length=8, color='k')

    ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax1.legend(prop={'size': 20},loc='best')#,shadow=True)

    ax1.set_ylim(ylim[0],ylim[1])
    ax1.set_xlim(xlim[0],xlim[1])
    ax1.set_xlabel(label[0], fontsize=25)#, fontdict=dict(weight='bold'))
    ax1.set_ylabel(label[1], fontsize=25)#, fontdict=dict(weight='bold'))
    fontsize=25
    for tick in ax1.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

    for tick in ax1.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        
    plt.savefig(pltname, bbox_inches = "tight")


# In[287]:


PrintFigures=False


# # NEGF

# ### parameter definition

# In[293]:


h=1.06e-34
q=1.602e-19
G_quantum=2*q**2/h #at zero temperature
print('G_quantum',G_quantum)
gamma=2.5 #eV
acc=0.142e-9
m=18
n=0
nano=1e-9

#N=num_layers_ch
#U=avg_field_ch

def get_radius(acc,m,n):
    return acc*(np.sqrt(3.)/(2*np.pi))*np.sqrt(m**2 + m*n + n**2)


R = get_radius(acc,m,n)
print('Radius of CNT (nm):', R/nano)

Eg_min = acc*gamma/R
print('bandgap, H/(eV):', Eg_min)

Us=0.0 #V, Source potential
Ud=0.0 #V, Drain potential

N=4
U=np.zeros(N)
print(U)
print('number of layers in the channel: ', N)

q = math.gcd(m, n)
Int=(m-n)/(3*q)
print('q',q)
print('Int',Int)
if(Int.is_integer()):
    d_R = 3*q
else:
    d_R = q
print('d_R',d_R)
Natoms_per_UC=int(4*(m**2+m*n+n**2)/d_R)
print('number of atoms per unit cell',Natoms_per_UC)

if(n==0): #zigzag
    M =  int(Natoms_per_UC/4)
elif(m==n): #armchair
    M = int(Natoms_per_UC/2)
else:
    M = np.abs(m-n)
print('number of total modes or atoms along the ring', M)
M_red = M
print('number of modes used to construct hamiltonian', M_red)


# ### Bandgap (From Theory)

# In[294]:


def get_bandgap_for_mode(p,acc,gamma,m,n):
    R=get_radius(acc,m,n)
    Eg_case=np.zeros(6,dtype=float)
    if(p>0):
        Eg_case[0] = (acc*gamma/R)*np.abs(3*p - (2*m+n))
        Eg_case[1] = (acc*gamma/R)*np.abs(3*p - (m+2*n))
        Eg_case[2] = (acc*gamma/R)*np.abs(3*p - (m-n))
        Eg_case[3] = (acc*gamma/R)*np.abs(-3*p + (2*m+n))
        Eg_case[4] = (acc*gamma/R)*np.abs(-3*p + (m+2*n))
        Eg_case[5] = (acc*gamma/R)*np.abs(3*p - (m-n))
    elif(p<0):
        Eg_case[0] = (acc*gamma/R)*np.abs(-3*p - (2*m+n))
        Eg_case[1] = (acc*gamma/R)*np.abs(-3*p - (m+2*n))
        Eg_case[2] = (acc*gamma/R)*np.abs(-3*p - (m-n))
        Eg_case[3] = (acc*gamma/R)*np.abs(3*p + (2*m+n))
        Eg_case[4] = (acc*gamma/R)*np.abs(3*p + (m+2*n))
        Eg_case[5] = (acc*gamma/R)*np.abs(-3*p - (m-n))
     
    Eg = np.amin(Eg_case)
    min_index = np.argmin(Eg_case)
    return (Eg,min_index+1)


m_arranged=np.arange(1,M+1,1)
m_arranged_rev=np.zeros(np.size(m_arranged),dtype=int)
counter=0
for i in reversed(m_arranged):
    m_arranged_rev[counter] = -i
    counter=counter+1
M_arr = np.concatenate((m_arranged_rev,m_arranged))
#M_arr = np.arange(-m+1,m,1)

print('M_arr',M_arr)

Eg_arr = np.zeros(M, dtype=float)
counter=0
for p in range(1,M+1):
    eg, case = get_bandgap_for_mode(p,acc,gamma,m,n)
    Eg_arr[p-1]=eg
    print('mode, bandgap', p , Eg_arr[p-1], case)
    #counter +=1

Eg_min = np.amin(Eg_arr)
print('minimum bandgap, acc*gamma/R',Eg_min, acc*gamma/R)

Eg_min_index_arr = np.where(Eg_arr==Eg_arr.min())
print('band with minimum bandgap:')
for j in Eg_min_index_arr:
    print(M_arr[M+j])
    
sorted_M_arr = [x for _,x in sorted(zip(Eg_arr,np.arange(1,M+1,1)))] 
#print('sorted_M_arr:',sorted_M_arr)
print('size of sorted_M_arr:',np.shape(sorted_M_arr))

M_red_arr = np.zeros(M_red,dtype=int)

for i in range(0,M_red):
    M_red_arr[i] = sorted_M_arr[i]
    
print('M_red_arr',M_red_arr)


# ### Band structure and First Brillouin Zone
# 

# In[295]:


def get_dispersion(m,gamma,nu,E_nu,kxa_arr):
    kyb=2*np.pi*nu/(2*m)
    for kxa in kxa_arr:
        E_nu.append(gamma*np.sqrt(1+(4*np.cos(kyb)*np.cos(kxa))+(4*np.cos(kyb)**2)))

k_pi=np.linspace(-np.pi,np.pi,2401) # range*pi    
k = k_pi/np.pi

nu_label=str(sorted_M_arr)
#print(nu_label)
E_k=[]
linewidth=[]
counter=0
lw=3


for p in sorted_M_arr:
    E_nu=[]
    get_dispersion(m,gamma,p,E_nu,k_pi)
    E_k.append(E_nu)
    E_nu_inv = []
    for i in range(0,len(k)):
        E_nu_inv.append(-1*E_nu[i])
    E_k.append(E_nu_inv)
    linewidth.append(lw)
    linewidth.append(lw)
    counter += 1
    if (counter%2 == 0):
        lw=4*np.exp(-0.2*counter)

print('size of E_k',np.shape(E_k))

if(PrintFigures):
    custom_plot_single(2*M*[k], 
                       E_k,
                       2*M*[None],
                       [-1,1],[-8,8],[r'k$_x$ a/$\pi$',r'Energy, E / (eV)'],'Dispersion_all.png',
                       2*M*['k'],
                       2*M*['solid'],
                       2*M*[None],
                       2*M*['none'],
                       linewidth,2*M*['none'])


# In[296]:


import operator
p1=sorted_M_arr[0]
p2=sorted_M_arr[1]
print('p2',p2)

Eg = gamma*acc/R
print('gamma',gamma,acc,R)
a=3*acc/2 
b=np.sqrt(3)*acc/2

def get_linear_approximate_dispersion(m,gamma,nu,E_nu,kxa_arr,kf):
    kyb=2*np.pi*nu/(2*m)
    factor=-1
    for kxa in kxa_arr:
        if(kxa > 0):
            factor =1
        E_nu.append((3*gamma/2)*np.sqrt((kxa*2/3-factor*kf[0]*acc)**2 + ((kyb*2/np.sqrt(3)-kf[1]*acc))**2))      

kf=[np.pi/a, np.pi/(3*b)]
E_lin_1 = []
get_linear_approximate_dispersion(m,gamma,p1,E_lin_1,k*np.pi,kf)
E_lin_1_minus=list(map(operator.mul, E_lin_1, [-1]*len(E_lin_1))) #multiply all elements in list by -1

kf=[0, 2*np.pi/(3*b)]
E_lin_2 = []
get_linear_approximate_dispersion(m,gamma,p2,E_lin_2,k*np.pi,kf)
E_lin_2_minus=list(map(operator.mul, E_lin_2, [-1]*len(E_lin_2))) #multiply all elements in list by -1

if(PrintFigures):
    custom_plot_single([k,k,k,k,k,k,k,k], 
                       [E_k[2*1-2],E_k[2*1-1],E_k[2*2-2],E_k[2*2-1],
                        E_lin_1,E_lin_1_minus,E_lin_2, E_lin_2_minus],
                       [r'$\nu$='+str(int(p1)),None, r'$\nu$='+str(int(p2)), None, r'$\nu$='+str(int(p1))+' (linear)',None, r'$\nu$='+str(int(p2))+' (linear)',None],
                       [-1,1],[-8,8],[r'k$_x$ a/$\pi$',r'Energy, E / (eV)'],'Dispersion.png',
                       ['k','k','b','b','r','r','g','g'],
                       ['solid','solid','solid','solid','dashed','dashed','dashed','dashed'],
                       8*[None],
                       8*['none'])
print('minimum band gap: (eV)', Eg)


# ### Construct hamiltonian H

# In[297]:


#The following 2 functions are to create off-diagonal elements of a hamiltonian
def get_zeros_and_ones(L,lst):
    for i in range(0,L):
        if(i%2==0):
            lst.append(0)
        else:
            lst.append(1)
            
def get_ones_and_zeros(L, lst): 
    for i in range(0,L):
        if(i%2==0):
            lst.append(1)
        else:
            lst.append(0)

lst_01 = [] 
lst_10 = [] 

get_zeros_and_ones(N-1, lst_01)
get_ones_and_zeros(N-1, lst_10)

def get_beta(gamma, M, J): #overlap coefficient in mode-space approximation
    return 2*gamma*np.cos(np.pi*J/M) #*np.exp(-1j*np.pi*J/M)


# In[298]:


H_mode=np.zeros((M_red,N,N),dtype=float)
print(np.shape(H_mode))
zero_mat=np.zeros((N,N),dtype=float)
#print(np.shape(zero_mat))
H=np.zeros((M_red*N,M_red*N),dtype=float)
print(np.shape(H))
beta = np.zeros(M_red,dtype=float)

for i in range(0,M_red):
    p=M_red_arr[i]
    beta[i] = get_beta(gamma,M,p)
    H_mode[i] =U*np.diag(np.ones(N)) + gamma*np.diag(lst_01,-1)  + beta[i]*np.diag(lst_10,-1) \
    + gamma*np.diag(lst_01,1) + beta[i]*np.diag(lst_10,1) 

H_lst=[]
for i in range (0, M_red):
    submat = []
    for j in range(0, M_red):
        if j==i:
            submat.append(H_mode[i])
        else:
            submat.append(zero_mat)          
    H_lst.append(submat)
H=np.block(H_lst)

if(PrintFigures):
    plt.figure(figsize = (16,16))
    sns.heatmap(H,cmap='coolwarm',annot=True)


# ### define E-space
# 

# In[299]:


Emin=-4.001 #eV
Emax=4.001 #eV
Epts=5000

E=np.linspace(Emin,Emax,Epts,dtype=float)
#print('E',E)

#Define fermi energies and potential due to bias
Ef_arr = np.linspace(-1, 1,100,dtype=float)
Cond_tot = np.zeros(np.size(Ef_arr),dtype=float)

kb=8.617333262e-5 #eV/K
Temp=298 #K
print(kb*Temp)

for f in range(0, np.size(Ef_arr)):
    Ef = Ef_arr[f]
    mu_s=Ef + Us
    mu_d=Ef + Ud 
    mu_c = Ef
    f_s=1./(1.+np.exp((E-mu_s)/(kb*Temp)))
    f_d=1./(1.+np.exp((E-mu_d)/(kb*Temp)))
    f_c=1./(1.+np.exp((E-mu_c)/(kb*Temp)))
    print(mu_s, mu_d, mu_c)
    
    dfdE_s = np.gradient(f_s)
    dfdE_d = np.gradient(f_d)
    dfdE_c = np.gradient(f_c)
    
    
    # In[300]:
    if(PrintFigures):
        custom_plot_single([E-Ef,E-Ef,E-Ef], 
                           [ f_s,f_d,f_c],
                           [r'source (%.2f V)'%(Us),r'drain (%.2f V)'%(Ud),r'channel (%.f K)'%(Temp)],
                           [-0.5,0.5],[0,1],[r'Energy, E / (eV)',r'Fermi function, f'],'fermi_functions.png',
                           ['k','r','g'],
                           ['solid','dashed','solid','dashed'],
                           ['o','^','o','^'],
                           4*['none'],
                           [2,2,2,2])
    
        custom_plot_single([E-Ef], 
                           [-dfdE_c ],
                           [r'channel (%.f K)'%(Temp)],
                           [-0.5,0.5],[0,0.005],[r'Energy, E / (eV)',r'Gradient of Fermi function, -$\frac{\partial{f}}{\partial{E}}$'],'gradient_of_fermi_functions.png',
                           ['k','r','k','r'],
                           ['solid','dashed','solid','dashed'],
                           ['o','^','o','^'],
                           4*['none'],
                           [2,2,2,2])
    
    
    # ### compute self-energy
    # 
    
    # In[301]:
    
    
    def get_analytical_retarted_surface_GF(E, Us,gamma,beta_j,Print):
        EmUs = E-Us
        Factor=EmUs**2 + gamma**2 - beta_j**2
        Sqrt=cmath.sqrt(Factor**2 - 4 * EmUs**2 * gamma**2)
        Denom = 2 * gamma**2 * EmUs
    
        Numer1= (Factor + Sqrt)
        Numer2= (Factor - Sqrt)
      
        if(Print):
            print('E,Us,gamma,beta_j',E,Us,gamma,beta_j)
            print('*',4 * EmUs**2 * gamma**2/Factor**2)
            print('Sqrt',Sqrt)
            print('Denom',Denom)
            print('Numer1',Numer1)
            print('Numer2',Numer2)
            print(np.imag(Numer1),np.imag(Numer2))
    
        zplus=1j*1e-5 #small imaginary number
    
        val1 = Numer1/Denom - zplus
        val2 = Numer2/Denom - zplus
        
        val = 0. + 1j*0.
        if(np.imag(val1) < 0.):
            val = val1 
        elif (np.imag(val2) < 0.):
            val = val2 
            
        if(Print): print(val)
    #     if(np.abs(np.imag(val)) < 1e-2):
    #         val = -zplus
        return val
    
    
    # In[302]:
    
    
    gR_surf_s1 = np.zeros((M_red,Epts),dtype=complex)
    gR_surf_d1 = np.zeros((M_red,Epts),dtype=complex)
    Ud=0
    print('Us,Ud',Us,Ud)
    
    for j in range(0,M_red):
        Print=False
        for e in range(0,Epts):
            #if(j==10):
            #    Print=True
            gR_surf_s1[j][e] = get_analytical_retarted_surface_GF(E[e],Us,gamma,beta[j],Print)
            gR_surf_d1[j][e] = get_analytical_retarted_surface_GF(E[e],Ud,gamma,beta[j],Print)
    
    
    # In[303]:
    
    
    if(PrintFigures):
        custom_plot_single([E,E], 
                           [ np.imag(gR_surf_s1[0]),np.imag(gR_surf_d1[0]) ],
                           [r'$g_1$, source ('+str(Us)+' V)',r'$g_1$, drain ('+str(Ud)+' V)'],
                           [-9,9],[-0.45,0.1],[r'Energy, E / (eV)',r'Surface Green\'s function, $g_1$'],'g1.png',
                           ['r','r','b','b'],
                           ['solid','dashed','solid','dashed'],
                           4*[None],#['o','^','o','^'],
                           4*['none'],
                           [2,3,2,3])
    print('minimum band gap/2 : (eV)', Eg/2)
    
    
    # ### compute Self-energy and Level-broadning
    
    # In[304]:
    
    
    Sigma_s_comp = np.zeros((M_red,Epts),dtype=complex)
    Sigma_d_comp = np.zeros((M_red,Epts),dtype=complex)
    Gamma_s_comp = np.zeros((M_red,Epts),dtype=complex)
    Gamma_d_comp = np.zeros((M_red,Epts),dtype=complex)
    
    for j in range(0,M_red):
        Sigma_s_comp[j] = gamma**2 * gR_surf_s1[j]
        Sigma_d_comp[j] = gamma**2 * gR_surf_d1[j]
            
        Gamma_s_comp[j] = 1j*(Sigma_s_comp[j] - np.conjugate(Sigma_s_comp[j]))
        Gamma_d_comp[j] = 1j*(Sigma_d_comp[j] - np.conjugate(Sigma_d_comp[j]))
    
    
    # In[305]:
    
    
    if(PrintFigures):
        custom_plot_single([E,E,E,E], 
                           [ np.real(Sigma_s_comp[0]), np.imag(Sigma_s_comp[0]), np.real(Sigma_d_comp[0]), np.imag(Sigma_d_comp[0])],
                           [r'$\Re(\Sigma_s^{1,1})$',r'$\Im(\Sigma_s^{1,1})$',r'$\Re(\Sigma_d^{N,N})$',r'$\Im(\Sigma_d^{N,N})$'],
                           [0,4],[-3,3],[r'Energy, E / (eV)',r'Component of self-energy, $\Sigma$'],'Sigma.png',
                           ['r','r','b','b'],
                           ['solid','dashed','solid','dashed'],
                           ['o','^','o','^'],
                           4*['none'],
                           [2,3,2,3])
    print('minimum band gap/2 : (eV)', Eg/2)
    
    
    # In[306]:
    
    
    if(PrintFigures):
        custom_plot_single([E,E], 
                           [ np.real(Gamma_s_comp[0]), np.real(Gamma_d_comp[0])],
                           [r'$\Gamma_s^{1,1}$',r'$\Gamma_d^{N,N}$'],
                           [-8,8],[0,5],[r'Energy, E / (eV)',r'Level broadning, $\Gamma$'],'Gamma.png',
                           ['r','b'],
                           ['solid','dashed'],
                           2*[None],
                           2*['none'],
                           [3,3])
    
    
    # In[307]:
    
    
    from numpy.linalg import inv
    zplus=1j*1e-12 #small imaginary number
    
    T=np.zeros(Epts,dtype=float)
    Cond=np.zeros(Epts,dtype=float)
    
    D=np.zeros(Epts,dtype=float)
    
    Cond_tot[f] = 0.
    I=0. #Current
    deltaE=E[1]-E[0]
    print(deltaE)
    
    
    # In[308]:
    
    
    for e in range(0,Epts):
        Sigma_s=np.zeros((M_red*N,M_red*N),dtype=complex)
        Sigma_d=np.zeros((M_red*N,M_red*N),dtype=complex)
        Gamma_s=np.zeros((M_red*N,M_red*N),dtype=complex)
        Gamma_d=np.zeros((M_red*N,M_red*N),dtype=complex)
        for j in range(0,M_red):
            Sigma_s[j*N][j*N] = Sigma_s_comp[j][e]
            Sigma_d[(j+1)*N-1][(j+1)*N-1] = Sigma_d_comp[j][e]
            Gamma_s[j*N][j*N] = Gamma_s_comp[j][e]
            Gamma_d[(j+1)*N-1][(j+1)*N-1] = Gamma_d_comp[j][e]
    #         print('\Sigma_s',Sigma_s[j*N][j*N], Sigma_s_comp[j][e])
    #         print('\Sigma_d',Sigma_d[(j+1)*N-1][(j+1)*N-1], Sigma_d_comp[j][e])
    #         print('\Gamma_s',Gamma_s[j*N][j*N], Gamma_s_comp[j][e])
    #         print('\Gamma_d',Gamma_d[(j+1)*N-1][(j+1)*N-1], Gamma_d_comp[j][e])  
    
        #or do the following (OK because Sigma matrix is diagonal)
        #Gamma_s = 1j*(Sigma_s - Sigma_s.conj().T) 
        #Gamma_d = 1j*(Sigma_d - Sigma_d.conj().T) 
        #sns.heatmap(np.real(Sigma_d+Sigma_s),cmap='coolwarm',annot=True)
        
        G_R=inv(((E[e]+zplus)*np.eye(M_red*N))- H - Sigma_s - Sigma_d)
        G_A = G_R.conj().T
        T[e]=np.real(np.trace(Gamma_s@G_R@Gamma_d@G_A))
        Cond[e]= -G_quantum*T[e]*dfdE_c[e]
        D[e]=np.real(np.trace(1j*(G_R - G_A)/(2*np.pi)))
        Cond_tot[f] += Cond[e]*deltaE
    
         
    #    I=I+(dE*IE*T[e]*(f1[k]-f2[k])) 
    
    
    # In[309]:
    
    
    print('Total conductance, Ef',Cond_tot[f], Ef)
    
    
    # In[310]:
    
    
    if(PrintFigures):
        custom_plot_single([E], 
                           [T/M_red],
                           [r'$T(E) (%d,%d)$'%(m,n)],
                           [-8,10],[0,1.1],[r'Energy, E / (eV)',r'Transmission probability, T'],'Transmission_Prob_%d_%d.png'%(m,n),
                           ['k'],
                           ['solid'],
                           1*[None],
                           1*['none'],
                           [3])
    
    
    # In[311]:
    
    
    print(Cond)
    
    
    # In[312]:
    
    
    if(PrintFigures):
        custom_plot_single([E], 
                           [Cond/G_quantum],
                           [r'$(%d,%d)$'%(m,n)],
                           [-2,2],[0,1e-6],[r'Energy, E / (eV)',r'Integrand of Conductance function, $G/(2q^2/h)$'],'Conductance_%d_%d.png'%(m,n),
                           ['k'],
                           ['solid'],
                           1*[None],
                           1*['none'],
                           [3])
    
    
    # In[313]:
    
    
    print(D)
    
    
    # # Comparison with Mintmire & White (Theory and First principles)
    
    # In[314]:
    
    
    nu0=np.round(2*m/3)
    print(nu0+M)
    #Mode_arr_Mintmire=np.linspace(-M,M,2*M+1)
    Mode_arr_Mintmire=np.linspace(-M,M,2*M+1)
    
    #Mode_arr_datta=np.linspace(nu0-M,nu0+M,2*M+1)
    print('Mode_arr_Mintmire',Mode_arr_Mintmire)
    #print('Mode_arr_datta',Mode_arr_datta)
    
    Lambda_Mintmire = 2*R/acc
    print('Lambda_Mintmire',Lambda_Mintmire)
    
    print('sorted_M_arr',sorted_M_arr)
    
    
    # In[315]:
    
    
    ##Francois' book
    D_th_lin = np.zeros(Epts, dtype=float)
    
    for e in range(0,Epts):
        const = acc*np.sqrt(3)/(np.pi**2*R*gamma)
        Sum=0.
        for nu in Mode_arr_Mintmire:   
            eps_m = np.abs(3*nu-m+n)*acc*gamma/(2*R) 
            if(np.abs(E[e]) > np.abs(eps_m)):
                Sum += np.abs(E[e])/cmath.sqrt(E[e]**2 - eps_m**2)
            elif(np.abs(E[e]) < np.abs(eps_m)):
                Sum +=0
        D_th_lin[e] = const*np.real(Sum)#*num_atoms
    
    dE_dk = np.zeros((M,np.size(k)),dtype=float)
    for p in range(1,M+1):
        dE_dk[p-1] = np.gradient(E_k[p*2-2])
        
    # D_th = np.zeros(Epts, dtype=float)   
    # for e in range(0,Epts):
    #     print('e',e)
    #     const = acc**2*np.sqrt(3)/(np.pi**2*R)
    #     Sum=0.+1j*0.
    #     for p in range(1,M+1):
    #         #print(p)
    #         k_sum=0.
    #         for kx in range(0,np.size(k)):
    #             convergence=np.abs((np.abs(E[e])-np.abs(E_k[p*2-2][kx]))/E[e])
    #             if(convergence < 1e-3): 
    #                 k_sum += 1./(np.abs(dE_dk[p-1][kx]))
    #                 if(np.abs(dE_dk[p-1][kx])==0):
    #                     print('mode, convergence, E, E(k)', sorted_M_arr[p-1], convergence, E[e],E_k[p*2-2][kx])
    #             else:
    #                 k_sum=0.
    #         Sum += k_sum*(k[1]-k[0])
    #     D_th[e] = const*np.real(Sum)
    
    
    ###Supriyo Datta's
    #
    # D_th_datta = np.zeros(Epts, dtype=float)
    
    # for e in range(0,Epts):
    #     a=3*acc/2.
    #     L=1*(3/4)*acc/17
    #     const = (2*L/(np.pi*a*gamma))
    #     Sum=0.
    #     for nu in Mode_arr_datta:
    #         Ek = (gamma*2*np.pi/np.sqrt(3))*(3*nu/(2*m)-1)+1j*1e-12
    #         Sum += E[e]/cmath.sqrt(E[e]**2 - Ek**2)
            
    #     D_th_datta[e] = const*Sum
    
    
    # In[316]:
    
    
    if(PrintFigures):
        custom_plot_single([Lambda_Mintmire*E/gamma,Lambda_Mintmire*E/gamma], 
                           [Lambda_Mintmire*gamma*D*2/M/N,Lambda_Mintmire*gamma*D_th_lin],
                           [r'NEGF [(%d,%d) CNT]'%(m,n), r'Mintmire & White (Theory)'],
                           [-5,5],[0,10],[r'Normalized energy, $\Lambda E$ / $\gamma$',r'Normalized Density of States, $\Lambda \gamma D(E)$'],'DOS_MintmireCompare_%d_%d.png'%(m,n),
                           ['k','r','b'],
                           ['solid','dashed','dotted'],
                           3*[None],
                           3*['none'],
                           3*[3])
    
    
    # In[317]:
    
    
    def get_bandgap_from_DOS(E, Dnorm):
        Epts = np.size(E)
        e_plus_traverse = True
        e_minus_traverse = True
        e_plus = 0.
        e_minus = 0.
        bandgap = 0.
        if(Dnorm[int(Epts/2.)] > 1e-3):
            bandgap = 0. 
        else:
            for e in range(int(Epts/2.), Epts):
                if(e_plus_traverse):
                    if(E[e] > 0.):
                        if(Dnorm[e] > 0.1):
                            e_plus = E[e]
                            e_plus_traverse = False
                            break
                    
            for e in np.arange(int(Epts/2.), 0,-1):
                if(e_minus_traverse):
                    if(E[e] < 0.):
                        if(Dnorm[e] > 0.1):
                            e_minus = E[e]
                            e_minus_traverse = False
                            break
                            print(E[e],D[e])
            bandgap = e_plus - e_minus
        return bandgap
    
    bandgap_negf = get_bandgap_from_DOS(E, Lambda_Mintmire*gamma*D*2/M/N )
    print('m, n, bandgap_negf, bandgap_theory: ', m, n, bandgap_negf, Eg_min)


fig, ax1 = plt.subplots(1,1)
fig.set_size_inches(10,10)
fig.patch.set_facecolor('white')
ax1.patch.set_facecolor('white')
ax1.plot(Ef_arr, Cond_tot, linewidth=3,  color='k', marker='o',fillstyle='none')
ax1.legend(prop={'size': 20},loc='best')
ax1.tick_params(which='minor', width=2, length=4, color='k')
ax1.tick_params(which='major', width=2, length=8, color='k')

ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator())
ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax1.legend(prop={'size': 20},loc='best')#,shadow=True)

ax1.set_yscale('log')
ax1.set_xlim(Ef_arr[0], Ef_arr[-1])
ax1.set_xlabel(r'Fermi level of the channel, Ef / (eV)', fontsize=25)#, fontdict=dict(weight='bold'))
ax1.set_ylabel(r'Conductance, $G$ / (S)', fontsize=25)#, fontdict=dict(weight='bold'))
fontsize=25
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)

for tick in ax1.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

plt.savefig('Total_conductance_vs_FermiLevel.png', bbox_inches = "tight")
