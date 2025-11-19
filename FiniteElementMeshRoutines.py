import numpy as np
import math
#from geopy import distance
import re


deg2km=111.132954
deg2rad=np.pi/180

def Distance(x0,y0,x1,y1):
#    d = distance.great_circle(y0,x0,y1,x1).km
    d=np.sqrt( 
        ( deg2km*np.cos(deg2rad*(y0+y1)/2 )*(x1-x0) )**2 + 
        ( deg2km*(y1-y0))**2  ) 
    return d


def loadWW3MeshCoords(fl):
    f=open(fl, 'r')
    header = f.readline() 
    header = f.readline() 
    header = f.readline() 
    header = f.readline() 
    header = f.readline() # number of nodes
    nn = re.findall(r'\d+', header)
    nn=int(nn[0])
    print(nn)
    xi=np.zeros(nn)
    yi=np.zeros(nn)
    k=0
    for i in range(nn):
        A = f.readline()
        #print(A)
        values = A.split(" ")
        #print(values)
        if len(values)>4:
            xi[k]=values[2]
            yi[k]=values[4]
        else:
            xi[k]=values[1]
            yi[k]=values[2]
#        print(str(xi[k])+" "+str(yi[k]))
        k=k+1
    print("number of nodes read: "+str(k))
    header = f.readline() 
    header = f.readline() 
    header = f.readline() # number of elements
    ne=int(header)
    print("ne=",str(ne))
    ei=np.zeros((ne,3), dtype=int)
    #print(ei)
    k=0
    for i in range(ne):
        A = f.readline()
        #print(A)
        values = A.split(" ")
        #print(values)
        #need to factor bnd nodes
        if len(values)>15:
            ei[k,0]=int(values[12])
            ei[k,1]=int(values[14])
            ei[k,2]=int(values[16])
            k=k+1
        elif len(values)>7:
            ei[k,0]=int(values[6])
            ei[k,1]=int(values[7])
            ei[k,2]=int(values[8])
            k=k+1
    print("number of elements read: "+str(k))
    return xi, yi, ei


def lengthscale(x, y, e):
    ne=e.shape[0]
    lengthscaleE=np.zeros(ne)
    for k in range(ne):
        if k % 10000 ==0:
            print(k)
        x1=x[e[k,0]-1]
        y1=y[e[k,0]-1]
        x2=x[e[k,1]-1]
        y2=y[e[k,1]-1]
        x3=x[e[k,2]-1]
        y3=y[e[k,2]-1]
        D3=Distance(x1,y1,x2,y2)
        D1=Distance(x2,y2,x3,y3)
        D2=Distance(x3,y3,x1,y1)
        lengthscaleE[k]=(D1+D2+D3)/3  
        # mean edge length
    return lengthscaleE

def ElementArea(x, y, e):
    ne=e.shape[0]
    AreaE=np.zeros(ne)
    for k in range(ne):
        if k % 10000 ==0:
            print(k)
        x1=x[e[k,0]-1];
        y1=y[e[k,0]-1];
        x2=x[e[k,1]-1];
        y2=y[e[k,1]-1];
        x3=x[e[k,2]-1];
        y3=y[e[k,2]-1];
        D3=Distance(x1,y1,x2,y2)
        D1=Distance(x2,y2,x3,y3)
        D2=Distance(x3,y3,x1,y1)
        S=(D1+D2+D3)/2
        A=np.sqrt(S * (S - D1) * (S - D2) * (S - D3))
        AreaE[k]=A
    return AreaE

def ComputeNodeLengthScale(x,y,e):
    nn=np.max(e[:])
    print("ComputeNodeLengthScale nn="+str(nn))
    areaE=ElementArea(x, y, e)
    LengthScaleE=lengthscale(x, y, e)
    
    ne=e.shape[0]
    areaT=np.zeros(nn)
    lsN=np.zeros(nn)
    for k in range(ne):
        for j in range(3):
            i=e[k,j]-1
            lsN[i]=lsN[i]+LengthScaleE[k]*areaE[k]
            areaT[i]=areaT[i]+areaE[k]
    for k in range(nn):
        lsN[k]=lsN[k]/areaT[k]
    return lsN

def ComputeNodeLengthScale0(lsE, areaE, e):
    nn=np.max(e[:])
    print("ComputeNodeLengthScale nn="+str(nn))
    ne=e.shape[0]
    lsT=np.zeros(nn)
    areaT=np.zeros(nn)
    lsN=np.zeros(nn)
    for k in range(ne):
        for j in range(3):
            i=e[k,j]-1
            lsT[i]=lsT[i]+lsE[k]*areaE[k]
            areaT[i]=areaT[i]+areaE[k]
    for k in range(nn):
        lsN[k]=lsT[k]/areaT[k]
    return lsN

