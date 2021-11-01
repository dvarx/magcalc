import numpy as np
import matplotlib.pyplot as plt
import os
from math import floor,sqrt,pi,sin,cos
from numpy.linalg import norm
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
    os.chdir(os.getcwd()+"\\python_comps")
    N=13  # no field samples in each direction
    filenames=["bxs_wocore.csv","bys_wocore.csv","bzs_wocore.csv"]

    system=magnetic_system(["bxs_wocore.csv","bys_wocore.csv","bzs_wocore.csv"],"comsol_csv",N)

    #calculate field
    currents=np.array([10,10,10,-10,-10,-10])           #opposing fields
    #currents=np.array([10,10,10,-10,-10,-10])          #additive fields
    Bsum=system.getBact((0,0,0)).dot(currents)

    system.plot_yz_view(currents)

    print("finished")