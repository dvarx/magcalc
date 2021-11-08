# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 13:15:09 2021

@author: dvarx
"""

import plotly.graph_objects as go
import numpy as np
import pickle
from numpy.linalg import norm

fptr=open("comsol/6a_add_xmag_xyz.pkl","rb")
(usig,usint,currents_app,Bs)=pickle.load(fptr)
fptr.close()

N=usig.shape[1]

zs_sweep=np.linspace(-2e-2,2e-2,N)
ys_sweep=np.linspace(-2e-2,2e-2,N)
xs_sweep=np.linspace(-2e-2,2e-2,N)

values=np.zeros((N,N,N))
X=np.zeros((N,N,N))
Y=np.zeros((N,N,N))
Z=np.zeros((N,N,N))

idxx=5
idxy=5
idxz=5

ref_sig=usig[:,idxx,idxy,idxz]/norm(usig[:,idxx,idxy,idxz])

for k in range(0,N):
    for j in range(0,N):
        for i in range(0,N):
            X[i,j,k]=xs_sweep[i]
            Y[i,j,k]=ys_sweep[j]
            Z[i,j,k]=zs_sweep[k]
            values[i,j,k]=1-(usig[:,i,j,k]/norm(usig[:,i,j,k])).dot(ref_sig)


fig = go.Figure(data=go.Isosurface(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    value=values.flatten(),
    surface_count=10,
    opacity=0.75,
    caps=dict(x_show=False, y_show=False)
    ))

fig.show()
