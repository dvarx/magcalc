import numpy as np
import matplotlib.pyplot as plt
import os
from math import floor,sqrt,pi,sin,cos
from numpy.linalg import norm
#from mayavi import mlab
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
    currents=np.array([10,10,10,10,10,10])          #opposing fields
    #currents=np.array([10,10,10,-10,-10,-10])        #additive fields
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
    #interpolator for a  magnetic field
    Bfield_int=getInterpolator(Bsum)
    #interpolator for a magnetic actuation matrix
    Bact_int=getInterpolator(Bs)

    #compute the null space of the actuation matrix in the center
    U,S,V = np.linalg.svd(Bact_int((0,0,0)))


    #search the edge of the null space for the optimal current
    Bact=Bact_int((0,0,0))
    ds=10e-3
    Bact_plusdx=Bact_int((ds,0,0))
    Bact_minusdx=Bact_int((-ds,0,0))
    Bact_plusdy=Bact_int((0,ds,0))
    Bact_minusdy=Bact_int((0,-ds,0))
    Bact_plusdz=Bact_int((0,0,ds))
    Bact_minusdz=Bact_int((0,0,-ds))
    #helper function for determning the gradient
    def getGradient(i):
            #compute minimum field growth in any direction
            dB_per_ds_x=min(norm(Bact_plusdx.dot(i))/ds,norm(Bact_minusdx.dot(i))/ds)
            dB_per_ds_y=min(norm(Bact_plusdy.dot(i))/ds,norm(Bact_minusdy.dot(i))/ds)
            dB_per_ds_z=min(norm(Bact_plusdz.dot(i))/ds,norm(Bact_minusdz.dot(i))/ds)
            dB_per_ds_min=min((dB_per_ds_x,dB_per_ds_y,dB_per_ds_z))
            return dB_per_ds_min
    #compute a reference solution for comparison
    iref=np.array([10,10,10,10,10,10])
    dB_per_ds_min_ref=getGradient(iref)
    #the allowable current space is null(Bact) ∩ ||i||2 < sqrt(dim(null(Bact)))*imax
    imax=sqrt(3)*10
    N=100
    phis=np.linspace(0,2*pi,N)
    thetas=np.linspace(0,pi/2,N)
    #look at the field B(0+dr) next to the focal point and determine the field norm there. for each current,
    #store the field norm in the worst case direction, e.g. the direction that generates the smallest norm(B)
    #indices of result array dB_per_ds_arr : dB_per_ds_arr[0] ~ phi coordinate , dB_per_ds_arr[1] ~ theta coordinate , dB_per_ds_arr[2] ~ dB_per_ds_min
    def getNullSpaceEdgeCurrent(theta,phi):
            alpha=cos(theta)*cos(phi)*imax
            beta=cos(theta)*sin(phi)*imax
            gamma=sin(theta)*imax
            return alpha*V[3,:]+beta*V[4,:]+gamma*V[5,:]
    dB_per_ds_arr=np.zeros((N*N,3))
    idx=0
    for phi in phis:
        for theta in thetas:
            dB_per_ds_arr[idx,0]=phi
            dB_per_ds_arr[idx,1]=theta
            #compute the current vector in the allowable current space
            i=getNullSpaceEdgeCurrent(phi,theta)
            #compute the field in the center (for debugging purposes)
            B=Bact.dot(i)
            print("Norm of current solution: %f"%(np.linalg.norm(B)))
            #compute minimum field growth in any direction
            dB_per_ds_x=min(norm(Bact_plusdx.dot(i))/ds,norm(Bact_minusdx.dot(i))/ds)
            dB_per_ds_y=min(norm(Bact_plusdy.dot(i))/ds,norm(Bact_minusdy.dot(i))/ds)
            dB_per_ds_z=min(norm(Bact_plusdz.dot(i))/ds,norm(Bact_minusdz.dot(i))/ds)
            dB_per_ds_min=min((dB_per_ds_x,dB_per_ds_y,dB_per_ds_z))
            dB_per_ds_arr[idx,2]=dB_per_ds_min
            idx+=1

    #find the current combination with the maximum dB_per_ds
    maxidx=np.argmax(dB_per_ds_arr[:,2])
    maxcurrents=getNullSpaceEdgeCurrent(dB_per_ds_arr[maxidx,0],dB_per_ds_arr[maxidx,1])
    maxBfield_int=getInterpolator(Bs.dot(maxcurrents))

    #plot the B field alongside the x axis
    xs=np.linspace(-2.5e-2,2.5e-2,60)
    ys=xs
    zs=xs
    magBsxaxis=[1e3*np.linalg.norm(maxBfield_int((x,0,0))) for x in xs]
    magBsyaxis=[1e3*np.linalg.norm(maxBfield_int((0,y,0))) for y in ys]
    magBszaxis=[1e3*np.linalg.norm(maxBfield_int((0,0,z))) for z in zs]

    plt.plot(1e3*xs,magBsxaxis)
    plt.plot(1e3*ys,magBsyaxis)
    plt.plot(1e3*zs,magBszaxis)

    plt.legend(("x axis","y axis","z axis"))
    plt.ylabel("|B| [mT]")
    plt.xlabel("Distance [mm]")
    plt.title("Currents: [%.2f,%.2f,%.2f,%.2f,%.2f,%.2f]"%(tuple([maxcurrents[i] for i in range(0,6)])))
    plt.grid()
    plt.show()

    #mlab.quiver3d(coordinates[:,0],coordinates[:,1],coordinates[:,2],Bsum[:,:,:,0].flatten(),Bsum[:,:,:,1].flatten(),Bsum[:,:,:,2].flatten())

    #mlab.show()

    print("finished")