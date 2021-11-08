# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 08:59:52 2021

@author: dvarx
"""
import numpy as np
from math import gcd,pi
import matplotlib.pyplot as plt

def lcm(numbers):
    lcm = 1
    for i in numbers:
        lcm = lcm*i//gcd(lcm, i)
    return lcm

f0=5e3      #the highest frequency in the system
Nlcm=lcm([7,8,9])
dT=1/f0/7
T1=9*dT
T2=8*dT
T3=7*dT
Trep=Nlcm*dT
f1=1/T1
f2=1/T2
f3=1/T3
Nsamples=Nlcm*100
ts=np.linspace(0,Trep,Nsamples)
dt=Trep/(Nsamples)
iampl=6
i1s=iampl*np.sin(2*pi*f1*ts)
i2s=iampl*np.sin(2*pi*f2*ts)
i3s=iampl*np.sin(2*pi*f3*ts)
plt.figure(1)
ax = plt.axes(projection='3d')
ax.plot3D(i1s,i2s,i3s)