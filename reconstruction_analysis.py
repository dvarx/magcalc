# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 16:44:37 2021

@author: dvarx
"""

import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
import pickle

fptr=open("comsol/6a_add_xmag_xyz.pkl","rb")
(usig,usint,currents_app,Bs)=pickle.load(fptr)
fptr.close()

#data extracted here defines the plot plane
#zy plane
usig=usig[:,5,:,:]
usint=usint[:,5,:,:]
axes=["y","z"]
#zy plane
#usig=usig[:,:,5,:]
#usint=usint[:,:,5,:]
#axes=["x","z"]

N=usig.shape[1]

zs_sweep=np.linspace(-2e-2,2e-2,N)
ys_sweep=np.linspace(-2e-2,2e-2,N)

Z=np.zeros((N,N))
Y=np.zeros((N,N))
VAL=np.zeros((N,N))

test_iy=8
test_iz=3
#we compare the other mpi signals with this reference signal
ref_sig=usig[:,test_iy,test_iz]/norm(usig[:,test_iy,test_iz])

for j in range(0,N):
    for i in range(0,N):
        Y[i,j]=ys_sweep[i]
        Z[i,j]=zs_sweep[j]
        VAL[i,j]=1-(usig[:,i,j]/norm(usig[:,i,j])).dot(ref_sig)
        
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})


ax.scatter(Y[test_iy,test_iz],Z[test_iy,test_iz],VAL[test_iy,test_iz],color="red")
ax.set_xlabel(axes[0])
ax.set_ylabel(axes[1])
ax.set_title(r"Cost Function $c(\vec{r})=(1-<\hat{s},s(\vec{r})>)^2$")

# Plot the surface.
surf = ax.plot_surface(Y, Z, VAL ,linewidth=0,alpha=0.75)

plt.show()