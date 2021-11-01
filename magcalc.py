import numpy as np
import matplotlib.pyplot as plt
import os
from math import floor,sqrt,pi,sin,cos,gcd,sinh,cosh
from numpy.linalg import norm
from scipy.interpolate import RegularGridInterpolator
from scipy.fft import fft

mu0=4*pi*1e-7
kb=1.38e-23
Tp=309
NA=6.022e23

def langevin(x):
    if x==0:
        return 0
    else:
        return cosh(x)/sinh(x)-1/x

def linidx2volidx(n):
    global N
    return (n%N,floor((n%N**2)/N),floor(n/N**2))

def volidx2linidx(nx,ny,nz):
    global N
    return nz*N**2+ny*N+nx

def get_magfield_data(filenames,N):
    Bs=np.zeros((N,N,N,3,6))    # array containing field data

    for coordidx in range(0,3):
        # open the file containing Bx, By or Bz values
        fptr=open(filenames[coordidx])
        linecounter=-1

        # read the magnetic flux components from the file containing Bx, By and Bz
        for line in fptr:
            linecounter+=1

            #ignore the first 5 header lines
            if linecounter<4:
                continue

            #third line contains coordinates
            if linecounter==4:
                split_lines=line.split("\"")
                coordinates=np.zeros((N**3,3))
                
                coordcounter=0
                for n in range(0,len(split_lines)):
                    #the odd line fragment contain the coordinates
                    if n%2==1:
                        #extract coordinate strings
                        coord_strs=split_lines[n][48:-1].split(",")
                        coordinates[coordcounter,:]=np.array([float(coord_strs[0]),float(coord_strs[1]),float(coord_strs[2])])
                        coordcounter+=1
                        if coordcounter>=N**3-1:
                            break

            #the lines contain the actual simulation data
            if linecounter>4:
                split_frags=line.split(",")
                for n in range(0,len(split_frags)):
                    if n<6:
                        pass
                    else:
                        linidx=n-6
                        nx=linidx2volidx(linidx)[0]
                        ny=linidx2volidx(linidx)[1]
                        nz=linidx2volidx(linidx)[2]
                        coilidx=linecounter-5
                        if(nz>=N):
                            break
                        Bs[nx,ny,nz,coordidx,coilidx]=float(split_frags[n])

        fptr.close()
    
    

    return (Bs,coordinates)

def lcm(numbers):
    lcm = 1
    for i in numbers:
        lcm = lcm*i//gcd(lcm, i)
    return lcm

class magnetic_system:
    def __init__(self,calibration_data,calibration_type,N):
        if(calibration_type=="comsol_csv"):
            self.Bs,self.coordinates=get_magfield_data(calibration_data,N)
            #extract the coordinate values of the grid
            self.xcoords=list(set(list(self.coordinates[:,0])))
            self.xcoords.sort()
            self.xmin=min(self.xcoords)
            self.xmax=max(self.xcoords)
            self.ycoords=list(set(list(self.coordinates[:,1])))
            self.ycoords.sort()
            self.ymin=min(self.xcoords)
            self.ymax=max(self.xcoords)
            self.zcoords=list(set(list(self.coordinates[:,2])))
            self.zcoords.sort()
            self.zmin=min(self.xcoords)
            self.zmax=max(self.xcoords)
            self.BactInterpolator=RegularGridInterpolator((self.xcoords,self.ycoords,self.zcoords),self.Bs)
        elif(calibration_type=="pickle"):
            raise Exception("Method not implemented")
        else:
            raise ValueError("calibration data need to be either comsol generated csv file or pickled calibration data")

    def getBact(self,pos):
        return self.BactInterpolator(pos)

    def plot_yz_view(self,currents):
        N=20
        zs=np.linspace(self.zmin,self.zmax,N)
        ys=np.linspace(self.zmin,self.zmax,N)
        vals=np.zeros((N,N))
        for i in range(0,N):
            for j in range(0,N):
                vals[i,j]=norm((self.BactInterpolator((0,ys[j],zs[i]))).dot(currents))
        plt.imshow(vals)
        plt.show()

if __name__=="__main__":
    # extract field data
    os.chdir(os.getcwd()+"\\comsol")
    N=13  # no field samples in each direction
    filenames=["\\comsol\\bxs_wocore.csv","\\comsol\\bys_wocore.csv","\\comsol\\bzs_wocore.csv"]

    system=magnetic_system(["bxs_wocore.csv","bys_wocore.csv","bzs_wocore.csv"],"comsol_csv",N)

    #calculate field
    #currents=np.array([10,10,10,-10,-10,-10])           #opposing fields
    #currents=np.array([10,10,10,-10,-10,0])          #additive fields
    #Bsum=system.getBact((0,0,0)).dot(currents)

    #system.plot_yz_view(currents)

    ###########################
    #mpi simulation
    ###########################

    #compute the currents for the six coils
    #-------------------------------------------------------------------------------------
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

    #define the concentration function & magnetization function
    #-------------------------------------------------------------------------------------
    Msat=6.5e3              #saturation magnetization [A/m]
    V=4*pi/3*(30e-9/2)**3   #particle volume [m^3]
    m=Msat*V                #maximum magnetic dipole moment of particle [Am^2]
    beta=mu0*m/(kb*Tp)
    #c=Msat/m                #particle concentration for unity iron concentration
    c=116*NA/(1e-3)         #166mol/l as reported by Panagiotopoulos1 et al
    #computes the magnetization M[A/m] dependent upon external magnetic field H[A/m] and iron volumetric concentration [m^-3]
    alpha=100  #unknown correction factor
    def M(c,H):
        global mu0,m,beta
        return c*m*langevin(alpha*beta*H)
    Bs=np.linspace(0,20e-3,100)
    Ms=np.zeros(Bs.shape)
    for i in range(0,len(Bs)):
        Ms[i]=M(c,Bs[i]/mu0)
    plt.figure()
    plt.plot(1e3*Bs,Ms)
    plt.xlabel("B [mT]")
    plt.show()

    #define the phantom concentration function
    #-------------------------------------------------------------------------------------

    #main simulation
    #-------------------------------------------------------------------------------------
    #particle positions to be simulated
    Nsweep=10
    zs_sweep=np.linspace(-2e-2,2e-2,Nsweep)
    ys_sweep=np.linspace(-2e-2,2e-2,Nsweep)
    usint=np.zeros((len(ts),Nsweep,Nsweep))
    usig=np.zeros((len(ts),Nsweep,Nsweep))
    Bs=np.zeros((len(ts),Nsweep,Nsweep))
    
    for sx in range(0,Nsweep):
        for sy in range(0,Nsweep):
            #compute the magnetic field at the origin
            zp=zs_sweep[sx]
            yp=ys_sweep[sy]
            print("Simulating signals with particle at position [%f,%f,%f]"%(0,yp,zp))
            for i in range(1,len(ts)):
                t=i*dt
                currents=np.array([i1s[i],i2s[i],i3s[i],i1s[i],i2s[i],i3s[i]])
                Bsum=system.getBact((0,yp,zp)).dot(currents)
        
                #compute the saturation factor (value of langevin function)
                Bs[i,sx,sy]=norm(Bsum)
                usint[i,sx,sy]=langevin(alpha*beta*norm(Bsum)/mu0)
                usig[i,sx,sy]=(usint[i,sx,sy]-usint[i-1,sx,sy])/dt


    #select signals to plot
    usig_plt=usig[:,5,5]
    usint_plt=usint[:,5,5]
    Bs_plt=Bs[:,5,5]

    fig, axs = plt.subplots(3,sharex=True)
    fig.suptitle('MPI Signals')
    axs[0].plot(ts*1e3,usint_plt)
    axs[0].set_ylabel(r"$\mathcal{L}(M)=\frac{M}{M_{sat}}$")
    axs[0].grid()
    axs[1].plot(ts*1e3,usig_plt)
    axs[1].grid()
    axs[1].set_ylabel(r"$\frac{d}{dt}\mathcal{L}(t)$")
    axs[2].plot(ts*1e3,Bs_plt)
    axs[2].grid()
    axs[2].set_ylabel(r"$|B(\vec{r_{P}},t)|$")
    plt.show()
    
    #compute and plot fft
    fig, axs = plt.subplots(2,sharex=True)
    #compute current ffts
    i1s_fft=fft(i1s)
    i2s_fft=fft(i2s)
    i3s_fft=fft(i3s)
    fig.suptitle("Harmonic Contents")
    Nfft=floor(Nsamples/2)
    axs[0].plot(abs(i1s_fft[:Nfft]+i2s_fft[:Nfft]+i3s_fft[:Nfft]))
    axs[0].set_ylabel("Current Signal")
    #compute
    usig_fft=fft(usig)
    axs[1].plot(abs(usig_fft[:Nfft]))
    axs[1].set_ylabel("Particle Signal")
    plt.show()

    #compute the magnetic field for the pickup coil (sensitivity function)

    print("finished")