
import numpy as np
import re
covmod="Sph"
K=2
variance=1.
deg2km=111.132954
deg2rad=np.pi/180
UseMean=True
MeanComp="Mean"#,"Median","Zero"

def InverseDistance(x,y,z,xi,yi,K):
    nx=len(x)
    d=np.zeros(nx)
    for j in range(nx):
        d[j] = Distance(x[j],y[j],xi,yi)
    w=1.0/(d**K)
    zi=np.dot( z, w )/np.sum( w )
    #print("zi")
    #print(zi)
    return zi

def GaussMarkov(x,y,z,xi,yi,LengthScale,K,Vobs,V,MeanTyp):
    nx=len(x)
    f = np.zeros(nx)
    C = np.zeros((nx,nx))
    for j in range(nx):
        for k in range(nx):
            d = Distance(x[j],y[j],x[k],y[k])
            C[j,k]=CovarianceDistance(d,covmod,LengthScale,V,2)
        C[j,j]=C[j,j] + Vobs
        d = Distance(x[j],y[j],xi,yi)
        f[j]=CovarianceDistance(d,covmod,LengthScale,V,2)
    detC = np.linalg.det(C)
    if np.isclose(detC, 0):
        zi=float("inf")
        print("Cov matrix is likely singular.")
    else:
        match MeanTyp:
            case "Zero": # simple kriging type
                mu=0.
            case "Mean":
                mu=np.mean(z) # regional mean
            case "Median":
                mu=np.median(z) # regional median
            case "Nearest": # nearest value
                mu=z[np.argmax(f)] # max correlate==closest point for process mean
#                d0 = np.zeros(nx)
#                for j in range(nx):
#                    d0[j]=Distance(x[j],y[j],xi,yi)
#                mu=z[np.argmin(d0)] # closest point for process mean
            case _:
                mu=0.
        w = np.linalg.solve(C, f)
        zi=mu+np.dot(w,z-np.array(mu))
    return zi

def GaussMarkovUnkMean(x,y,z,xi,yi,LengthScale,K,Vobs,V):
    nx=len(x)
    f = np.zeros(nx+1)
    C = np.zeros((nx+1,nx+1))
    for j in range(nx):
        for k in range(nx):
            d = Distance(x[j],y[j],x[k],y[k])
            C[j,k]=CovarianceDistance(d,covmod,LengthScale,V,2)
        C[j,j]=C[j,j] + Vobs
        d = Distance(x[j],y[j],xi,yi)
        f[j]=CovarianceDistance(d,covmod,LengthScale,V,2)
        C[nx,j]=1
        C[j,nx]=1
    C[nx,nx]=0
    f[nx]=1
    detC = np.linalg.det(C)
    if np.isclose(detC, 0):
        zi=float("inf")
        print("Cov matrix is likely singular.")
    else:
        w = np.linalg.solve(C, f)
        wp = np.zeros(nx)
        for j in range(nx):
            wp[j]=w[j]
        zi= np.dot(wp,z)
    return zi

#from geopy import distance
def Distance(x0,y0,x1,y1):
#    d = distance.great_circle(y0,x0,y1,x1).km
    d=np.sqrt( 
        ( deg2km*np.cos(deg2rad*(y0+y1)/2 )*(x1-x0) )**2 + 
        (deg2km*(y1-y0))**2  ) 
    return d

def CovarianceDistance(d,model,ls,v,k):
    if model == "Exp":
        c = v*np.exp(- (d/ls )**k )
    if model == "Sph":
        c = v*( 1 - 1.5*(d/ls) + .5*((d/ls)**3)   )
        if d > ls:
            c=0
    return c
        
