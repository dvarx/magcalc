# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 16:44:37 2021

@author: dvarx
"""

import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt

usig=np.load("sim_data/usig.npy")
usint=np.load("sim_data/usint.npy")

N=usig.shape[1]

zs_sweep=np.linspace(-2e-2,2e-2,N)
ys_sweep=np.linspace(-2e-2,2e-2,N)

Z=np.zeros((N,N))
Y=np.zeros((N,N))
VAL=np.zeros((N,N))

test_ix=7
test_iy=7
#we compare the other mpi signals with this reference signal
ref_sig=usig[:,test_ix,test_iy]/norm(usig[:,test_ix,test_iy])

for i in range(0,N):
    for j in range(0,N):
        Y[i,j]=ys_sweep[j]
        Z[i,j]=zs_sweep[i]
        VAL[i,j]=1-(usig[:,i,j]/norm(usig[:,i,j])).dot(ref_sig)
        
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})


ax.scatter(Y[test_ix,test_iy],Z[test_ix,test_iy],VAL[test_ix,test_iy],color="red")
ax.set_xlabel("y")
ax.set_ylabel("z")
ax.set_title(r"Cost Function $c(\vec{r})=1-<\hat{s},s(\vec{r})>$")

# Plot the surface.
surf = ax.plot_surface(Z, Y, VAL ,linewidth=0,alpha=0.75)

plt.show()