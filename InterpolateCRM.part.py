import os
import argparse
import time
import numpy as np
import netCDF4 as nc
import csv

#import jigsawpy
import sys
from scipy.interpolate import RegularGridInterpolator
import re
#import ComputeMeshToMeshInterpWeights as mshint
import math
#from geopy import distance
import itertools
#from geopy.distance import haversine
#from pykrige.ok import OrdinaryKriging

import FiniteElementMeshRoutines as FE
import GaussMarkov as GM

coords_1 = (34.0522, -118.2437)  # Los Angeles
coords_2 = (40.7128, -74.0060)   # New York City

#Using Haversine
#distance_haversine = distance.haversine(coords_1, coords_2).km
#dist=distance.great_circle(coords_1, coords_2).km
#print(f"Haversine distance: {dist:.2f} km")
#maxval=100.


#[0->7] variables(dimensions): |S1 crs(), float64 lat(lat), float64 lon(lon), float32 z(lat, lon)
#[8->9]variables(dimensions): float64 x(x), float64 y(y), float32 z(y, x)
#[10]variables(dimensions): float64 lon(num_lon), float64 lat(num_lat), int16 bed_elevation(num_row, num_col)
"""
flnms=[
    "crm_vol1_2023.nc",
    "crm_vol2_2023.nc",
    "crm_vol3_2023.nc",
    "crm_vol4_2023.nc",
    "crm_vol5_2023.nc",
    "crm_vol7_2024.nc",
    "crm_vol9_2023.nc",
    "crm_vol10_2023.nc", 
    "crm_vol6_2023.nc",
    "crm_vol8_2023.nc",
    "RTopo_2_0_4_GEBCO_v2023_60sec_pixel.CRMformat.nc"
    ]
"""
flnms=[
    "crm_vol1_2023.nc.S250m.nc",
    "crm_vol2_2023.nc.S250m.nc",
    "crm_vol3_2023.nc.S250m.nc",
    "crm_vol4_2023.nc.S250m.nc",
    "crm_vol5_2023.nc.S250m.nc",
    "crm_vol7_2024.nc.S250m.nc",
    "crm_vol9_2023.nc.S250m.nc",
    "crm_vol10_2023.nc.S250m.nc", 
    "crm_vol6_2023.nc.S250m.nc",
    "crm_vol8_2023.nc.S250m.nc",
    "crm_southak.nc",# dx=700m, fillvalue= -2147483648, not used...
    "RTopo_2_0_4_GEBCO_v2023_60sec_pixel.CRMformat.nc",
#    "crm_socal_1as_vers2.nc.S250m.nc",# mostly redundent with 6 but some extra regions, need to synthesize
    ]
 

lambdaLL=.025
NpointsMax=1000
xlist=[]
ylist=[]
n=0
nxl=np.zeros(len(flnms))
nyl=np.zeros(len(flnms))
xmax=np.zeros(len(flnms))
xmin=np.zeros(len(flnms))
ymax=np.zeros(len(flnms))
ymin=np.zeros(len(flnms))
for fl in flnms:
#    print(f"'{fl}' has a length of {len(fl)}")
    data = nc.Dataset(fl,"r")
    x=np.array(data["lon"][:])
    x=x%360
    y=np.array(data["lat"][:])
    xlist.append(list(x))
    ylist.append(list(y))
    nxl[n]=len(x)
    nyl[n]=len(y)
    #setup quick lookup table for file exlusion
    xmax[n]=np.max(x)+lambdaLL*2
    xmin[n]=np.min(x)-lambdaLL*2
    ymax[n]=np.max(y)+lambdaLL*2
    ymin[n]=np.min(y)-lambdaLL*2
    n=n+1
"""
xlist=[]
ylist=[]
n=0
nxl=np.zeros(len(flnms))
nyl=np.zeros(len(flnms))
for fl in flnms:
#    print(f"'{fl}' has a length of {len(fl)}")
    data = nc.Dataset(fl,"r")
    x=np.array(data["lon"][:])
    x=x%360
    y=np.array(data["lat"][:])
    xlist.append(list(x))
    ylist.append(list(y))
    nxl[n]=len(x)
    nyl[n]=len(y)
    n=n+1
"""
def ZeroPadIntStr(N,K):
    ZPNs=str(N).zfill(K)
    return ZPNs

mshnm="RWPStest"
mesh="meshes/"+mshnm+".msh"
OutDir=mshnm+".files/"

xi, yi, ei = FE.loadWW3MeshCoords(mesh)

lsE=FE.lengthscale(xi, yi, ei)
areaE=FE.ElementArea(xi, yi, ei)
lsN=FE.ComputeNodeLengthScale(lsE, areaE, ei)

print(np.min(lsN))
print(np.mean(lsN))
print(np.max(lsN))
N=int(sys.argv[1])
Nzp=ZeroPadIntStr(N,5)
print("N= "+str(N)+" Nzp=",Nzp)
fln=OutDir+'NodeList.'+str(N)+'.txt'
floutGM=OutDir+'GM.'+Nzp+'.txt'
floutGMUnk=OutDir+'GMUnk.'+Nzp+'.txt'
floutID=OutDir+'InvDist.'+Nzp+'.txt'
floutNpts=OutDir+'NDataPoints.'+Nzp+'.txt'
floutClosest=OutDir+'ClosestValue.'+Nzp+'.txt'
#read in nodes to interpolate to
f=open(fln,"r")
header = f.readline()
header = f.readline() # number of nodes
nn0=int(header)
nn = re.findall(r'\d+', header)
nn=int(nn[0])
print("nn0 ="+str(nn0))
#make local node set to interpolate to
xil=np.zeros(nn)
yil=np.zeros(nn)
LSl=np.zeros(nn)
for k in range(nn):
    j = int(f.readline())
    xil[k]=xi[j]
    yil[k]=yi[j]
    LSl[k]=lsN[j]
    print(str(k)+" " +str(j)+" "+str(LSl[k]))

nn=len(xil)
nf=len(flnms)

#nn=103
ziD=np.zeros(nn)
ziSK=np.zeros(nn)
ziOK=np.zeros(nn)
ziClosest=np.zeros(nn)
ziOKerr=np.zeros(nn)
NumPoints=np.zeros(nn)

t0 = time.time()
for n in range(nn):
    xp=xil[n]%360
    yp=yil[n]
    LSp=2.0*LSl[n]#LocalLengthScale for interpolation
    xs=[]
    ys=[]
    zs=[]
    si=0
    for j in range(nf):
        if all([ xp < xmax[j],xp > xmin[j],yp < ymax[j],yp > ymin[j]]):
            x=np.array(xlist[j][:])
            y=np.array(ylist[j][:] )
            fl=flnms[j]
            jx=np.array(np.where( np.abs(xp-x) < lambdaLL ))
            jx=list(itertools.chain.from_iterable(jx))
            jy=np.array(np.where( np.abs(yp-y) < lambdaLL ))
            jy=list(itertools.chain.from_iterable(jy))
            data = nc.Dataset(fl,"r")
            for kx in jx:
                for ky in jy:
                        if kx and ky: 
                            z=data["z"][ky,kx]
                            if not z.mask:
                                zd=float(z.data)
                                xs.append(x[kx])
                                ys.append(y[ky])
                                zs.append(zd)
    ms=np.column_stack((xs,ys))                               
    msu, indx = np.unique(ms, axis=0, return_index=True)
    n0=ms.shape[0]
    nu0=msu.shape[0]
    if nu0 < n0:
        print("points not unique for node: "+str(n))
        print(ms.shape)
        print(msu.shape)
        print(indx)
        print(indx.shape)
        xs=np.array(xs)
        ys=np.array(ys)
        zs=np.array(zs)
        xs=xs[indx]
        ys=ys[indx]
        zs=zs[indx]
    if n  % 100 == 0:
        print("interpolating to:"+str(yp)+":"+str(xp))
        t1 = time.time()
        time_per_iter = (t1 - t0) / (n + 1)
        time_remaining = (nn - n) * time_per_iter / 60
        print("Progress: "+str(n+1)+" of "+ str(nn))
        print(f"Time remaining: {time_remaining:.2f} minutes")
        print(f"Average time per node: {60*time_per_iter:.2f} seconds")
        print("node: "+str(n)+" uses "+str(len(xs))+" data points")                        
    
    Npoints=len(xs)
    NumPoints[n]=Npoints
    
    #limit total number of points in interpolation if something went wrong
    if Npoints > NpointsMax:
        print("node "+str(n)+" has "+str(Npoints)+" data points. taking nearest: "+str(NpointsMax))
        xsp=np.zeros(NpointsMax)
        ysp=xsp
        zsp=xsp
        D=GM.Distance(xs,ys,xp,yp)
        for k in range(NpointsMax):
            ki=np.argmin(D)
            xsp[k]=xs[ki]
            ysp[k]=ys[ki]
            zsp[k]=zs[ki]
            D[ki]=float('inf')
        xs=xsp
        ys=ysp
        zs=zsp
    
    ziD[n]=GM.InverseDistance(xs,ys,zs,xp,yp,2)
    stderr=np.mean(np.abs(zs))/100.
    #stderr=np.std(zs)/10.
    if stderr < 1:
        stderr = 1.
    ziOK[n]=GM.GaussMarkovUnkMean(xs, ys, zs, xp, yp,LSp, 2, stderr)
    ziSK[n]=GM.GaussMarkov(xs, ys, zs, xp, yp, LSp, 2, stderr)
    D=GM.Distance(xs,ys,xp,yp)
    jc=np.argmin(D)
    ziClosest[n]=zs[jc]
    
np.savetxt(floutID,  ziD, fmt='%.6f', delimiter='\n')
np.savetxt(floutGM,  ziSK, fmt='%.6f', delimiter='\n')
np.savetxt(floutGMUnk,  ziOK, fmt='%.6f', delimiter='\n')
np.savetxt(floutClosest,  ziClosest, fmt='%.6f', delimiter='\n')
np.savetxt(floutNpts, NumPoints, fmt='%d', delimiter='\n')
    
    
