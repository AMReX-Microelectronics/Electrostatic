#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
                       markevery=[45,30,50,40,56,72,63,95],
                       show_legend=True,
                       plt_outside=False):
    
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
    if(show_legend):
        ax1.legend(prop={'size': 20},loc='best')
    ax1.tick_params(which='minor', width=2, length=4, color='k')
    ax1.tick_params(which='major', width=2, length=8, color='k')

    ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    if(show_legend):
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
    if(plt_outside==False):
        plt.savefig(pltname, bbox_inches = "tight")
    return plt,pltname


# In[2]:


PrintFigures=False
PointCharge=True


# ### parameter definition

# In[3]:


h=1.06e-34
q=1.602e-19
G_quantum=2*q**2/h #at zero temperature
print('G_quantum',G_quantum)
print('R_quantum',1/G_quantum)
gamma=2.5 #eV
acc=0.142e-9
m=17
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

N=192
print('number of layers in the channel: ', N)

q_factor = math.gcd(m, n)
Int=(m-n)/(3*q_factor)
print('q_factor',q_factor)
print('Int',Int)
if(Int.is_integer()):
    d_R = 3*q_factor
else:
    d_R = q_factor
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
M_red = 1
print('number of modes used to construct hamiltonian', M_red)

nring_per_unitcell = 4 #for a zigzag nanotube
print('nring_per_unitcell',nring_per_unitcell)

avg_distance_per_layer = 3*acc/nring_per_unitcell
print('avg_distance_per_layer',avg_distance_per_layer)
print('Total length', N*avg_distance_per_layer)


# ### Point Charge

# In[4]:


Y_offset = -(N/2)*avg_distance_per_layer #nm
#print('Y_offset',Y_offset)
axial_loc = np.empty(N,dtype=float)

def get_axial_loc(r,ACC):
    y=0.
    if(r % 2 == 0):
        y = ACC
    elif (r%2 == 1):
        y = ACC/2
    return y
 
axial_loc[0]=acc/2.

for l in range(1,N):
    rID_in_unitcell = l%nring_per_unitcell   
#    print(rID_in_unitcell)
    axial_loc[l] = axial_loc[l-1] + get_axial_loc(rID_in_unitcell,acc)
    
for l in range(0,N):
    axial_loc[l] += Y_offset  
    
#print(axial_loc)

U=np.zeros(N)
pt_charge_loc=np.array([0.,1*nano])

Q=q #C
print('Q',q)
eps_r =1 
eps_0 = 8.85418782e-12 #F/m or C/(V.m)
if (PointCharge):
    for l in range(0, N):
        r = np.sqrt((axial_loc[l]-pt_charge_loc[0])**2 + pt_charge_loc[1]**2)
        U[l] = -1./(4.*np.pi*eps_0*eps_r) * Q/r
print('potential on layers:',U)
Us = U[0]
Ud = U[-1]

Ec0 = Eg_min/2
Ev0 = -Eg_min/2

print('Ec0',Ec0)
print('Ev0',Ev0)


# ### Bandgap (From Theory)

# In[5]:


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

M_red_arr = []
#np.zeros(M_red,dtype=int)

M_red_arr.append(sorted_M_arr[0])
i = 1
counter = 1
while(counter < M_red):
    if(Eg_arr[sorted_M_arr[i]-1] != Eg_arr[sorted_M_arr[i-1]-1]):
        M_red_arr.append(sorted_M_arr[i])
        counter +=1
    i += 1
print('M_red_arr',M_red_arr)
M_red_arr[0]=11
print('M_red_arr',M_red_arr)


# ### Band structure and First Brillouin Zone
# 

# In[6]:


def get_dispersion(m,gamma,nu,E_nu,kxa_arr):
    kyb=2*np.pi*nu/(2*m)
    for kxa in kxa_arr:
        E_nu.append(gamma*np.sqrt(1+(4*np.cos(kyb)*np.cos(kxa))+(4*np.cos(kyb)**2)))

k_pi=np.linspace(-np.pi,np.pi,2401) # range*pi    
k = k_pi/np.pi

nu_label=str(sorted_M_arr)
#print(nu_label)
E_k=[]
linecolor=[]
linewidth=[]
counter=0
lw=3
lc='r'

for p in sorted_M_arr:
    E_nu=[]
    get_dispersion(m,gamma,p,E_nu,k_pi)
    E_k.append(E_nu)
    E_nu_inv = []
    for i in range(0,len(k)):
        E_nu_inv.append(-1*E_nu[i])
    E_k.append(E_nu_inv)
    if(counter > 1):
        print('counter,p',counter,p)
        lc='k'  
    linecolor.append(lc)
    linecolor.append(lc)
    
    linewidth.append(lw)
    linewidth.append(lw)
    counter += 1
    if (counter%2 == 0):
        lw=4*np.exp(-0.2*counter)

print('size of E_k',np.shape(E_k))

if(PrintFigures):
    plt,pltname = custom_plot_single(2*M*[k], 
                       E_k,
                       2*M*[None],
                       [-1,1],[-8,8],[r'k$_x$ a/$\pi$',r'Energy, E / (eV)'],'Dispersion_all.png',
                       linecolor,
                       2*M*['solid'],
                       2*M*[None],
                       2*M*['none'],
                       linewidth,2*M*['none'],
                       show_legend=False,
                       plt_outside=True)
    plt.axhline(y = 0., color = 'b', linestyle = 'dashed',label='Fermi')
    plt.text(-0.7, 0.1, 'Fermi Level', fontsize=22,color='b')
    plt.savefig(pltname, bbox_inches = "tight")


# In[7]:


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
    plt,pltname=custom_plot_single([k,k,k,k,k,k,k,k], 
                       [E_k[2*1-2],E_k[2*1-1],E_k[2*2-2],E_k[2*2-1],
                        E_lin_1,E_lin_1_minus,E_lin_2, E_lin_2_minus],
                       [r'$\nu$='+str(int(p1)),None, r'$\nu$='+str(int(p2)), None, r'$\nu$='+str(int(p1))+' (linear)',None, r'$\nu$='+str(int(p2))+' (linear)',None],
                       [-1,1],[-8,8],[r'k$_x$ a/$\pi$',r'Energy, E / (eV)'],'Dispersion.png',
                       ['k','k','r','r','k','k','r','r'],
                       ['solid','solid','solid','solid','dashed','dashed','dashed','dashed'],
                       8*[None],
                       8*['none'],
                       plt_outside=True)
    
    plt.axhline(y = 0., color = 'b', linestyle = 'dashed',label='Fermi')
    plt.text(-0.7, 0.1, 'Fermi Level', fontsize=22,color='b')
    plt.savefig(pltname, bbox_inches = "tight")
    
print('minimum band gap: (eV)', Eg)


# ### Construct hamiltonian H

# In[8]:


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


# In[9]:


print(M_red_arr[0])
print(get_beta(gamma,M,11))


# In[10]:


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
#    sns.heatmap(H,cmap='coolwarm',annot=True)


# ### define E-space
# 

# In[11]:


Emin=-2.001 #eV
Emax=2.001 #eV
Epts=1000

E=np.linspace(Emin,Emax,Epts,dtype=float)

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
    print('mu_c',mu_c)
    f_s=1./(1.+np.exp((E-mu_s)/(kb*Temp)))
    f_d=1./(1.+np.exp((E-mu_d)/(kb*Temp)))
    f_c=1./(1.+np.exp((E-mu_c)/(kb*Temp)))
    print(mu_s, mu_d, mu_c)
    print('Us', Us)
    print('Ud',Ud)
    dfdE_s = np.gradient(f_s)
    dfdE_d = np.gradient(f_d)
    dfdE_c = np.gradient(f_c)
    
    
    # In[12]:
    
    
    if(PrintFigures):
        plt,pltname = custom_plot_single([axial_loc/nano,axial_loc/nano], 
                           [Ev0+U, Ec0+U],
                           [r'$Ev_0$-eVg', r'$Ec_0$-eVg'],
                           [axial_loc[0]/nano,axial_loc[-1]/nano],[-2,1],[r'z / (nm)',r'E / (eV)'],'Ec0_Ev0_minus_U.png',
                           ['b','r'],linestyle=['solid','solid'],
                           show_legend=True,
                           plt_outside=True)
        plt.axhline(y = Ef, color = 'k', linestyle = 'dashed',label='Fermi')
        plt.text(0, Ef, 'Fermi Level', fontsize=22,color='k')
    #    plt.text(-0.1, (Ev0), 'U=-eVg = 0.2 eV; Vg = -0.2 eV', fontsize=22,color='b')
    
    #    plt.axhline(y = Ev0, color = 'b', linestyle = 'dashed',label='Fermi')
    #    plt.text(-1, (Ev0), 'Ev0-U(z=10nm)', fontsize=22,color='b')
    
    #    plt.axhline(y = Ec0, color = 'r', linestyle = 'dashed',label='Fermi')
    #    plt.text(-1, (Ec0), 'Ec0-U(z=10nm)', fontsize=22,color='r')
    
    #    plt.axhline(y = Ev0-U[-1], color = 'b', linestyle = 'dashed',label='Fermi')
    #    plt.text(-1, (Ev0-U[-1]), 'Ev0-U(z=10nm)', fontsize=22,color='b')
    
    #     plt.axhline(y = Ec0-U[-1], color = 'r', linestyle = 'dashed',label='Fermi')
    #     plt.text(-9, (Ec0-U[-1]), 'Ec0-U(z=10nm)', fontsize=22,color='r')
    
    #     plt.axhline(y = np.min(Ev0-U), color = 'g', linestyle = 'dashed',label='Fermi')
    #     plt.text(-9,  np.min(Ev0-U), 'min(Ev0-U)', fontsize=22,color='g')
    
        plt.savefig(pltname, bbox_inches = "tight")
    
    
    # In[13]:
    
    
    if(PrintFigures):
        custom_plot_single([E-Ef,E-Ef,E-Ef], 
                           [ f_s,f_d,f_c],
                           [r'source (%.2f V)'%(Us),r'drain (%.2f V)'%(Ud),r'channel (%.f K)'%(Temp)],
                           [-0.5,0.5],[0,1],[r'Energy, (E - $E_f$) / (eV)',r'Fermi function, f'],'fermi_functions.png',
                           ['k','r','g'],
                           ['solid','dashed','solid','dashed'],
                           ['o','^','o','^'],
                           4*['none'],
                           [2,2,2,2])
    
        custom_plot_single([E-Ef], 
                           [-dfdE_c ],
                           [r'channel (%.f K)'%(Temp)],
                           [-0.5,0.5],[0,0.1],[r'Energy, (E - $E_f$) / (eV)',r'Gradient of Fermi function, -$\frac{\partial{f}}{\partial{E}}$'],'gradient_of_fermi_functions.png',
                           ['k','r','k','r'],
                           ['solid','dashed','solid','dashed'],
                           ['o','^','o','^'],
                           4*['none'],
                           [2,2,2,2])
    
    
    # ### compute self-energy
    # 
    
    # In[14]:
    
    
    def get_analytical_retarted_surface_GF(E, U,gamma,beta_j,Print):
        zplus = 1j*1e-14
    #    zplus=0
        EmU=E+zplus-U
    
        Factor=EmU**2 + gamma**2 - beta_j**2
        Sqrt=cmath.sqrt(Factor**2 - 4 * EmU**2 * gamma**2)
        Denom = 2 * gamma**2 * EmU
    
        Numer1= (Factor + Sqrt)
        Numer2= (Factor - Sqrt)
      
        if(Print):
            print('E,U,gamma,beta_j',E,U,gamma,beta_j)
            print('*',4 * EmUs**2 * gamma**2/Factor**2)
            print('Sqrt',Sqrt)
            print('Denom',Denom)
            print('Numer1',Numer1)
            print('Numer2',Numer2)
            print(np.imag(Numer1),np.imag(Numer2))
    
        zplus=1j*1e-8 #small imaginary number
    
        val1 = Numer1/Denom #- zplus
        val2 = Numer2/Denom #- zplus
    #     val = 0. + 1j*0.
    #     if(np.imag(val1) < 0.):
    #         val = val1 
    #     elif (np.imag(val2) < 0.):
    #         val = val2 
            
        if(Print): print(val)
        return val1
    
    
    # In[15]:
    
    
    print(Us,Ud)
    
    
    # In[16]:
    
    
    gR_surf_s1 = np.zeros((M_red,Epts),dtype=complex)
    gR_surf_d1 = np.zeros((M_red,Epts),dtype=complex)
    #Ud=0
    #print('Us,Ud',Us,Ud)
    
    for j in range(0,M_red):
        Print=False
        for e in range(0,Epts):
            #if(j==10):
            #    Print=True
            gR_surf_s1[j][e] = get_analytical_retarted_surface_GF(E[e],Us,gamma,beta[j],Print)
            gR_surf_d1[j][e] = get_analytical_retarted_surface_GF(E[e],Ud,gamma,beta[j],Print)
    
    
    # In[17]:
    
    
    if(PrintFigures):
        plt,pltname=custom_plot_single([ np.imag(gR_surf_s1[0]),np.imag(gR_surf_d1[0]) ],
                           [E,E], 
                           [r'$g_1$, source ('+str(Us)+' V)',r'$g_1$, drain ('+str(Ud)+' V)'],
                           [-0.45,0.1],[-2,2],[r'Surface Green\'s function, $g_1$',r'Energy, E / (eV)'],'g1.png',
                           ['r','r','b','b'],
                           ['solid','dashed','solid','dashed'],
                           4*[None],#['o','^','o','^'],
                           4*['none'],
                           [2,3,2,3], plt_outside=True)
        plt.axhline(y = Ev0+U[-1], color = 'b', linestyle = 'dashed',label='Fermi')
        plt.text(0.1, (Ev0+U[-1]), 'Ev0-eVg(z=10nm)', fontsize=22,color='b')
    
        plt.axhline(y = Ec0+U[-1], color = 'r', linestyle = 'dashed',label='Fermi')
        plt.text(0.1, (Ec0+U[-1]), 'Ec0-eVg(z=10nm)', fontsize=22,color='r')
    
        plt.savefig(pltname, bbox_inches = "tight")
    print('minimum band gap/2 : (eV)', Eg/2)
    
    
    # ### compute Self-energy and Level-broadning
    
    # In[18]:
    
    
    Sigma_s_comp = np.zeros((M_red,Epts),dtype=complex)
    Sigma_d_comp = np.zeros((M_red,Epts),dtype=complex)
    Gamma_s_comp = np.zeros((M_red,Epts),dtype=complex)
    Gamma_d_comp = np.zeros((M_red,Epts),dtype=complex)
    
    for j in range(0,M_red):
        Sigma_s_comp[j] = gamma**2 * gR_surf_s1[j]
        Sigma_d_comp[j] = gamma**2 * gR_surf_d1[j]
            
        Gamma_s_comp[j] = 1j*(Sigma_s_comp[j] - np.conjugate(Sigma_s_comp[j]))
        Gamma_d_comp[j] = 1j*(Sigma_d_comp[j] - np.conjugate(Sigma_d_comp[j]))
    
    
    # In[19]:
    
    
    if(PrintFigures):
    #     plt,pltname=custom_plot_single([ np.real(Sigma_s_comp[0]), np.imag(Sigma_s_comp[0]), np.real(Sigma_d_comp[0]), np.imag(Sigma_d_comp[0])],
    #                        [E,E,E,E], 
    #                        [r'$\Re(\Sigma_s^{1,1})$',r'$\Im(\Sigma_s^{1,1})$',r'$\Re(\Sigma_d^{N,N})$',r'$\Im(\Sigma_d^{N,N})$'],
    #                        [-3,3],[-2,2],[r'Component of self-energy, $\Sigma$',r'Energy, E / (eV)'],'Sigma.png',
    #                        ['r','r','b','b'],
    #                        ['solid','dashed','solid','dashed'],
    #                        ['o','^','o','^'],
    #                        4*['none'],
    #                        [2,3,2,3], plt_outside=True)
        plt,pltname=custom_plot_single([np.imag(Sigma_s_comp[0]),  np.imag(Sigma_d_comp[0])],
                           [E,E,E,E], 
                           [r'$\Im(\Sigma_s^{1,1})$',r'$\Im(\Sigma_d^{N,N})$'],
                           [-3,3],[-2,2],[r'Component of self-energy, $\Sigma$',r'Energy, E / (eV)'],'Sigma.png',
                           ['r','r','b','b'],
                           ['solid','dashed','solid','dashed'],
                           ['o','^','o','^'],
                           4*['none'],
                           [2,3,2,3], plt_outside=True)
        plt.axhline(y = Ev0+U[-1], color = 'b', linestyle = 'dashed',label='Fermi')
        plt.text(0.1, (Ev0+U[-1]), 'Ev0-eVg(z=10nm)', fontsize=22,color='b')
    
        plt.axhline(y = Ec0+U[-1], color = 'r', linestyle = 'dashed',label='Fermi')
        plt.text(0.1, (Ec0+U[-1]), 'Ec0-eVg(z=10nm)', fontsize=22,color='r')
    
        plt.savefig(pltname, bbox_inches = "tight")
    
    print('minimum band gap/2 : (eV)', Eg/2)
    
    
    # In[20]:
    
    
    if(PrintFigures):
        plt,pltname=custom_plot_single([np.real(Gamma_s_comp[0]), np.real(Gamma_d_comp[0])],
                           [E,E],
                           [r'$\Gamma_s^{1,1}$',r'$\Gamma_d^{N,N}$'],
                           [-0.1,5],[-2,2],[r'Level broadning, $\Gamma$',r'Energy, E / (eV)'],'Gamma.png',
                           ['r','b'],
                           ['solid','dashed'],
                           2*[None],
                           2*['none'],
                           [3,3], plt_outside=True)
        plt.axhline(y = Ev0+U[-1], color = 'b', linestyle = 'dashed',label='Fermi')
        plt.text(0.1, (Ev0+U[-1]), 'Ev0-eVg(z=10nm)', fontsize=22,color='b')
    
        plt.axhline(y = Ec0+U[-1], color = 'r', linestyle = 'dashed',label='Fermi')
        plt.text(0.1, (Ec0+U[-1]), 'Ec0-eVg(z=10nm)', fontsize=22,color='r')
    
    #    plt.axhline(y = np.min(Ev0-U), color = 'g', linestyle = 'dashed',label='Fermi')
    #    plt.text(2.1,  np.min(Ev0-U), 'min(Ev0-U)', fontsize=22,color='g')
    #    plt.text(1.25, 1.5, r'T=Trace[$\Gamma_s G^R \Gamma_d G^A$]', fontsize=22,color='k')
        plt.savefig(pltname, bbox_inches = "tight")
    
    
    # In[21]:
    
    
    from numpy.linalg import inv
    zplus=1j*1e-12 #small imaginary number
    
    T=np.zeros(Epts,dtype=float)
    Cond=np.zeros(Epts,dtype=float)
    
    D=np.zeros(Epts,dtype=float)
    A=np.zeros((Epts,N),dtype=float)
   
    Cond_tot[f] = 0.
    I=0. #Current
    deltaE=E[1]-E[0]
    print(deltaE)
    
    
    # In[22]:
    
    
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
            #H[j*N][j*N]=0
            #H[(j+1)*N-1][(j+1)*N-1]=0
    #         print('\Sigma_s',Sigma_s[j*N][j*N], Sigma_s_comp[j][e])
    #         print('\Sigma_d',Sigma_d[(j+1)*N-1][(j+1)*N-1], Sigma_d_comp[j][e])
    #         print('\Gamma_s',Gamma_s[j*N][j*N], Gamma_s_comp[j][e])
    #         print('\Gamma_d',Gamma_d[(j+1)*N-1][(j+1)*N-1], Gamma_d_comp[j][e])  
    
        #or do the following (OK because Sigma matrix is diagonal)
        #Gamma_s = 1j*(Sigma_s - Sigma_s.conj().T) 
        #Gamma_d = 1j*(Sigma_d - Sigma_d.conj().T) 
        #sns.heatmap(np.real(Sigma_d+Sigma_s),cmap='coolwarm',annot=True)
    #    if(PrintFigures):
    #        if(e == 0):
    #            plt.figure(figsize = (16,16))
    #            #sns.heatmap(np.real(Sigma_s+Sigma_d),cmap='coolwarm',annot=True)
    #            sns.heatmap(H,cmap='coolwarm',annot=True)
     
        G_R=inv(((E[e]+zplus)*np.eye(M_red*N))- H - Sigma_s - Sigma_d)
        G_A = G_R.conj().T
        T[e]=np.real(np.trace(Gamma_s@G_A@Gamma_d@G_R))
        Cond[e]= -G_quantum*T[e]*dfdE_c[e]
        Amat=1j*(G_R - G_A)
        D[e]=np.real(np.trace(Amat)/(2*np.pi))
        A[e] = np.real(np.diag(Amat))
        #print(A[0]) 
        Cond_tot[f] += Cond[e]*deltaE
    
         
    #    I=I+(dE*IE*T[e]*(f1[k]-f2[k])) 
    
    
    # In[23]:
    
    
    if(PrintFigures):
        plt,pltname=custom_plot_single([T], 
                           [E],
                           [r'$T(E) (%d,%d)$'%(m,n)],
                           [0,1.1],[-2,1],[r'Transmission, T(E)', r'Energy, E / (eV)', 'T'],'Transmission_%d_%d.png'%(m,n),
                           ['k'],
                           ['solid'],
                           1*[None],
                           1*['none'],
                           [3],show_legend=False,plt_outside=True)
        plt.axhline(y = Ev0+U[-1], color = 'b', linestyle = 'dashed',label='Fermi')
        plt.text(1.1, (Ev0+U[-1]), 'Ev0-eVg', fontsize=22,color='b')
    
        plt.axhline(y = Ec0+U[-1], color = 'r', linestyle = 'dashed',label='Fermi')
        plt.text(1.1, (Ec0+U[-1]), 'Ec0-eVg', fontsize=22,color='r')
    
        plt.axhline(y = np.min(Ev0+U), color = 'g', linestyle = 'dashed',label='Fermi')
        #plt.text(1.1,  np.min(Ev0+U), 'min(Ev0-U)', fontsize=22,color='g')
        #plt.text(0.5, 0.8, r'T=Trace[$\Gamma_s G^R \Gamma_d G^A$]', fontsize=22,color='k')
    
        plt.savefig(pltname, bbox_inches = "tight")
        
    
    
    # In[24]:
    
    
    #print(Cond)
    
    
    # In[25]:
    
    
    if(PrintFigures):
        plt,pltname=custom_plot_single([Cond/G_quantum], 
                           [E],
                           [r'$(%d,%d)$'%(m,n)],
                           [0,0.01],[-2,1],[r'Conductance function, $G(E)/(2q^2/h)$',r'Energy, E / (eV)'] ,'Conductance_%d_%d.png'%(m,n),
                           ['k'],
                           ['solid'],
                           1*[None],
                           1*['none'],
                           [3],show_legend=False,plt_outside=True)
        plt.axhline(y = Ef, color = 'k', linestyle = 'dashed',label='Fermi')
        plt.text(0.01, Ef, 'Fermi Level', fontsize=22,color='k')
        
        plt.axhline(y = Ev0+U[-1], color = 'b', linestyle = 'dashed',label='Fermi')
        plt.text(0.01, (Ev0+U[-1]), 'Ev0-U(z=10nm)', fontsize=22,color='b')
    
        plt.axhline(y = Ec0+U[-1], color = 'r', linestyle = 'dashed',label='Fermi')
        plt.text(0.01, (Ec0+U[-1]), 'Ec0-U(z=10nm)', fontsize=22,color='r')
    
        plt.axhline(y = np.min(Ev0+U), color = 'g', linestyle = 'dashed',label='Fermi')
        plt.text(0.01,  np.min(Ev0+U), 'min(Ev0-U)', fontsize=22,color='g')
        plt.text(0.004, 0.8, r'$G(E)=-\frac{2q^2}{h} T \frac{\partial f}{\partial E}$', fontsize=30,color='k')
        plt.text(0.004, 0.5, r'$G_{tot}=\int G(E) dE$', fontsize=30,color='k')
    
        plt.savefig(pltname, bbox_inches = "tight")
    
    
    # In[26]:
    
    
    #print(D)
    
    
    # In[27]:
    
    
    print('Total conductance:',Cond_tot)
    
    
    # In[28]:
    
    
    from matplotlib import ticker, cm
    import matplotlib.colors as colors
    def contourSpanVsTime(contourdata,ContourMin,ContourMax,xarr,yarr,fig, ax,cbarTrue,cbarlabel,
                         x_lst, y_lst, color_lst):
        print(len(xarr),len(yarr))
        print(xarr[0],xarr[-1])
        X, Y = np.meshgrid(len(xarr),len(yarr))
        cmap = plt.get_cmap('gist_heat')
        #levels = MaxNLocator(nbins=20).tick_values(0.0001, ContourMax)
        #levels=levels.astype(float)
    
        levels_exp = np.arange(np.floor(np.log10(1e-12)-1),np.ceil(np.log10(ContourMax)+1))
        levels = np.power(10, levels_exp)
        print(levels)
        extent=[xarr[0],xarr[-1],yarr[0],yarr[-1]]
    
        cpl=ax.contourf(contourdata,extent=extent,levels=levels,
                        norm=colors.LogNorm(vmin=1e-12, vmax=ContourMax),
                        cmap=cmap)
        for p in range(0, len(x_lst)):
            ax.plot(x_lst[p], y_lst[p], color_lst[p],linewidth=3,linestyle='solid')
    #    ax.clabel(cpl, inline=True, fontsize=15,fmt='%.1e',colors='r')
        if (cbarTrue==True):
            cbar = fig.colorbar(cpl,ax=ax,format='%.2e')
            #cbar.minorticks_on()
            cbar.ax.set_title(cbarlabel,fontsize = 25,weight="bold")
    
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax.tick_params(which='minor', length=3, color='black',width=1)
        ax.tick_params(which='major', length=5, color='black',width=1)
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        plt.yticks(fontsize=20)
        plt.xticks(fontsize=20)
    
        return plt
    
    from matplotlib import rc, rcParams
    rc('text', usetex=True)
    
    fig, ax = plt.subplots(1,1, constrained_layout=True)
    fig.set_size_inches(10,10)
    
    xe=N
    ye=Epts
    myplt1=contourSpanVsTime(A,A.min(),A.max(),axial_loc/nano,E,fig, ax, cbarTrue=True,cbarlabel=r'A / (eV$^{-1}$)',
                             x_lst=[axial_loc/nano,axial_loc/nano], 
                             y_lst=[Ev0+U, Ec0+U],
                             color_lst=['b','r'])
    ax.set_xlabel(r'z / (nm)',fontsize=25)
    ax.set_ylabel(r'E / (eV)',fontsize=25)
    plt.savefig('SpectralFunction.png', bbox_inches = "tight")
    
    
    # # Comparison with Mintmire & White (Theory and First principles)
    
    # In[29]:
    
    
    nu0=np.round(2*m/3)
    print(nu0+M)
    #Mode_arr_Mintmire=np.linspace(-M,M,2*M+1)
    Mode_arr_Mintmire=np.linspace(-M,M,2*M+1)
    
    #Mode_arr_datta=np.linspace(nu0-M,nu0+M,2*M+1)
#    print('Mode_arr_Mintmire',Mode_arr_Mintmire)
    #print('Mode_arr_datta',Mode_arr_datta)
    
    Lambda_Mintmire = 2*R/acc
#    print('Lambda_Mintmire',Lambda_Mintmire)
    
#    print('sorted_M_arr',sorted_M_arr)
    
    
    # In[30]:
    
    
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
        
    if(PrintFigures):
        custom_plot_single([Lambda_Mintmire*E/gamma,Lambda_Mintmire*E/gamma], 
                           [Lambda_Mintmire*gamma*D*2/M/N,Lambda_Mintmire*gamma*D_th_lin],
                           [r'NEGF [(%d,%d) CNT]'%(m,n), r'Mintmire & White (Theory)'],
                           [-5,5],[0,10],[r'Normalized energy, $\Lambda E$ / $\gamma$',r'Normalized Density of States, $\Lambda \gamma D(E)$'],'DOS_Normalized_%d_%d.png'%(m,n),
                           ['k','r','b'],
                           ['solid','dashed','dotted'],
                           3*[None],
                           3*['none'],
                           3*[3])
    
    
    if(PrintFigures):
        plt,pltname = custom_plot_single([Lambda_Mintmire*gamma*D*2/M/N],
                                         [E], 
                           [r'NEGF [(%d,%d) CNT]'%(m,n)],
                           [0,1],[-2,2],[r'Normalized Density of States, $\Lambda \gamma D(E)$',r'Energy, $E$ / (eV)'],'DOS_%d_%d.png'%(m,n),
                           ['k','r','b'],
                           ['solid','dashed','dotted'],
                           3*[None],
                           3*['none'],
                           3*[3], plt_outside=True)
        plt.text(0.3, 1, 'point charge= %de'%(Q/q), fontsize=22,color='k')
        plt.savefig(pltname, bbox_inches = "tight")
    
    

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

np.savetxt('Conductance_vs_FermiLevel_CNT_%d_%d_Temp_%d.dat'%(m,n,Temp),
        np.c_[Ef_arr,
              Cond_tot],
        header='"Ef_arr", "Gtot"',
        comments='', fmt='%.14e', delimiter=' ',newline='\n')
