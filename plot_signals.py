# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 10:53:18 2021

@author: dvarx
"""
import matplotlib.pyplot as plt
import numpy as np
import pickle
from math import gcd,floor
from numpy.linalg import norm
from scipy.fft import fft
from math import pi

fptr=open("comsol/6a_add_xmag.pkl","rb")
(usig,usint,currents_app,Bsmat,ts)=pickle.load(fptr)
fptr.close()

mu0=4*pi*1e-7
kB=1.38*1e-23
T=298

#select signals to plot
usig=usig[:,5,5]
usint=usint[:,5,5]
Bs=np.zeros(Bsmat.shape[1])
for i in range(0,Bsmat.shape[1]):
    Bs[i]=norm(Bsmat[:,i,5,5])
i1s=currents_app[0,:]
i2s=currents_app[1,:]
i3s=currents_app[2,:]

#scaling factors for induced voltage
#particle / concentration dependent scaling factors
mparticle=10e-6                     #particle mass [kg]
nparticles=6.2e17                   #particles per kg [#/kg]
Nparticles=nparticles*mparticle     #no of particles [#]
Dparticle=50e-9                     #diameter of iron core of SPIO
m=0.6/mu0*pi/6*Dparticle**3         #magnetic moment of particle (V_particle*Hsat)
#Hint, the upper value evalates to around 3.12e-17 [Am^2] whereas Goodwill and Conolly predict 3.1e-17 [Am^2] for 50nm particles
M=m*Nparticles                      #magnetic dipole moment of particles [Am^2]
Hx=177e-3                         #sensitivity of pickup coil in x direction

usig=-mu0*Hx*M*usig

fig, axs = plt.subplots(3,sharex=True)
fig.suptitle('MPI Signals')
axs[0].plot(ts*1e3,usint)
axs[0].set_ylabel(r"$\mathcal{L}(M)=\frac{M}{M_{sat}}$")
axs[0].grid()
axs[1].plot(ts*1e3,usig)
axs[1].grid()
axs[1].set_ylabel(r"$v_{pickup}$")
axs[2].plot(ts*1e3,Bs)
axs[2].grid()
axs[2].set_ylabel(r"$|B(\vec{r_{P}},t)|$")
axs[2].set_xlabel("Time t  [ms]")
plt.show()

#compute and plot fft
fig, axs = plt.subplots(2,sharex=True)
#compute current ffts
i1s_fft=fft(i1s)
i2s_fft=fft(i2s)
i3s_fft=fft(i3s)
fig.suptitle("Harmonic Contents")
Nsamples=ts.shape[0]
Nfft=floor(Nsamples/2)
axs[0].plot(abs(i1s_fft[:Nfft]+i2s_fft[:Nfft]+i3s_fft[:Nfft]))
axs[0].set_ylabel("Current Signal")
#compute
usig_fft=fft(usig)
axs[1].plot(abs(usig_fft[:Nfft]))
axs[1].set_ylabel("Particle Signal")
plt.show()

usig_rms=np.sqrt(1/len(usig)*usig.dot(usig))
print("RMS of inducd voltage signal: %f"%(usig_rms))