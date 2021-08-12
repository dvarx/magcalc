import numpy as np
import matplotlib.pyplot as plt
import os
from math import floor,sqrt,pi,sin,cos
from numpy.linalg import norm
from mayavi import mlab
from scipy.interpolate import RegularGridInterpolator

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
                        Bs[nx,ny,nz,coordidx,coilidx]=float(split_frags[n])

        fptr.close()
    
    return (Bs,coordinates)

if __name__=="__main__":
    # extract field data
    os.chdir("C:\\Users\\dvarx\\src\\magcalc\\python_comps")
    N=13                        # no field samples in each direction
    filenames=["bxs.csv","bys.csv","bzs.csv"]
    (Bs,coordinates)=get_magfield_data(filenames,N)

    #calculate field
    #currents=np.array([10,10,10,10,10,10])          #opposing fields
    currents=np.array([10,10,10,-10,-10,-10])        #additive fields
    Bsum=Bs.dot(currents)

    #extract the coordinate values of the grid
    xcoords=list(set(list(coordinates[:,0])))
    xcoords.sort()
    ycoords=list(set(list(coordinates[:,1])))
    ycoords.sort()
    zcoords=list(set(list(coordinates[:,2])))
    zcoords.sort()
    #define interpolator function for magnetic fields
    def getInterpolator(B):
        return RegularGridInterpolator((xcoords,ycoords,zcoords),B)

    field_int=getInterpolator(Bsum)

    mlab.quiver3d(coordinates[:,0],coordinates[:,1],coordinates[:,2],Bsum[:,:,:,0].flatten(),Bsum[:,:,:,1].flatten(),Bsum[:,:,:,2].flatten())

    mlab.show()

    print("finished")