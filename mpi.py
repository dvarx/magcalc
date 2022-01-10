# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 10:13:57 2022

@author: dvarx
"""

from math import pi
import numpy as np
from math import floor,sqrt,pi,sin,cos,gcd,sinh,cosh
import matplotlib.pyplot as plt

mu0=4*pi*1e-7
kb=1.38e-23
Tp=309
NA=6.022e23

def langevin(x):
    if x==0:
        return 0
    else:
        return cosh(x)/sinh(x)-1/x

#define the concentration function & magnetization function
#-------------------------------------------------------------------------------------
V=4*pi/3*(30e-9/2)**3               #particle volume [m^3]
Dparticle=50e-9
m=0.6/mu0*pi/6*Dparticle**3         #magnetic moment of particle (V_particle*Bsat)
beta=mu0*m/(kb*Tp)
cm=6.2e17               #specific particle concenration (#particles per kg)
msample=0.1e-3          #mass of sample [kg]


#computes the magnetization M[A/m] dependent upon external magnetic field H[A/m] and iron volumetric concentration [m^-3]
def mdipole(cm,msample,H):
    global mu0,m,beta
    return cm*msample*m*langevin(beta*H)

#plot magnetization curve of single dipole
#-------------------------------------------------------------------------------------
Bs=np.linspace(0,20e-3,100)
mdipoles=np.zeros(Bs.shape)
for i in range(0,len(Bs)):
    mdipoles[i]=mdipole(cm,msample,Bs[i]/mu0)
plt.figure(1)
plt.plot(1e3*Bs,mdipoles)
plt.xlabel("B [mT]")
plt.ylabel("Magnetic Moment [A/m]")
plt.show()

#compute response to sinusoidal excitation signal
#-------------------------------------------------------------------------------------
f0=10e3
T0=1/f0
Nsim=20
Tsim=Nsim*T0
ts=np.linspace(0,Tsim,Nsim*1000)
Bs=5e-3*np.sin(2*pi*f0*ts)
mdipoles=np.zeros(len(Bs))
for i in range(0,len(ts)):
    mdipoles[i]=mdipole(cm,msample,Bs[i]/mu0)
plt.figure(2)
plt.plot(ts,mdipoles)

#
#
dt=ts[1]-ts[0]
Hx=1e-3/mu0
vs=np.zeros(len(Bs))
vs[0]=0
for i in range(1,len(Bs)):
    vs[i]=mu0*Hx*(mdipoles[i]-mdipoles[i-1])/dt
plt.figure(3)
plt.plot(ts,vs)