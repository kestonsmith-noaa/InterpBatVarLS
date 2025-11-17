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
Zmin=-11000. #deepest legit ocean depth value in case of mask fail 
# nescesarry for data sets "crm_vol6_2023.nc.S250m.nc" and "crm_vol8_2023.nc.S250m.nc"
Zmax=0. # maximum value to include in interpolation - zero out land values or ignore land 
Dmin=10./1000.# minimum distance between observation points in km, prevent singularity
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

#Nzp=ZeroPadIntStr(N,5)
#print("N= "+str(N)+" Nzp=",Nzp)

fln=OutDir+'NodeList.'+str(N)+'.txt'
floutID=OutDir+'InvDist.'+str(N)+'.txt'
floutGMM=OutDir+'GMM.'+str(N)+'.txt'
floutGMN=OutDir+'GMN.'+str(N)+'.txt'
floutGM0=OutDir+'GM0.'+str(N)+'.txt'
floutGMU=OutDir+'GMU.'+str(N)+'.txt'
floutClosest=OutDir+'ClosestValue.'+str(N)+'.txt'
floutNpts=OutDir+'NDataPoints.'+str(N)+'.txt'
floutLLS=OutDir+'LengthScale.'+str(N)+'.txt'

#read in nodes to interpolate to
f=open(fln,"r")
header = f.readline()
nn =int( f.readline())

#make local node set to interpolate to
xil=np.zeros(nn)
yil=np.zeros(nn)
LSl=np.zeros(nn)
for k in range(nn):
    j = int(f.readline())
    xil[k]=xi[j]
    yil[k]=yi[j]
    LSl[k]=lsN[j]

nn=len(xil)
nf=len(flnms)

NumPoints=np.zeros(nn)
ziID=np.zeros(nn)
ziGMU=np.zeros(nn)
ziGMN=np.zeros(nn)
ziGMM=np.zeros(nn)
ziGM0=np.zeros(nn)
ziClosest=np.zeros(nn)
LocalLengthScale=np.zeros(nn)

t0 = time.time()
for n in range(nn):
    xp=xil[n]%360
    yp=yil[n]
    LSp=LSl[n]#LocalLengthScale for interpolation
    LocalLengthScale[n]=LSp
    xs=[]
    ys=[]
    zs=[]
    si=0
    for j in range(nf):
        if all([ xp < xmax[j],xp > xmin[j],yp < ymax[j],yp > ymin[j]]):
            x=np.array(xlist[j][:])
            y=np.array(ylist[j][:])
            lambdaLL0=lambdaLL
            fl=flnms[j]
            "RTopo_2_0_4_GEBCO_v2023_60sec_pixel.CRMformat.nc"
            if j==nf-1: # sparser dataset
                lambdaLL0=5.*lambdaLL # pull more points from sparser global dataset
            jx=np.array(np.where( np.abs(xp-x) < lambdaLL0 ))
            jx=list(itertools.chain.from_iterable(jx))
            jy=np.array(np.where( np.abs(yp-y) < lambdaLL0 ))
            jy=list(itertools.chain.from_iterable(jy))
            data = nc.Dataset(fl,"r")
            for kx in jx:
                for ky in jy:
                        if kx and ky: # not empty 
                            z=data["z"][ky,kx]
                            if not z.mask:
                                IncludePoint=True
                                if len(xs)>0:#BEGIN -check for near duplicate points - causes singularity in GM
                                    d=np.zeros(len(xs))
                                    for k in range(len(xs)):
                                        d[k]=GM.Distance(xs[k],ys[k],x[kx],y[ky])
                                    if np.min(d)<Dmin:
                                        print("near duplicate point at distance: "+str(np.min(d))+" for node "+str(n))
                                        IncludePoint=False #END -check for near duplicate points
                                zd=float(z.data)
                                if zd < Zmin:#double check mask fail and bad fill value
                                        IncludePoint=False
                                if IncludePoint: # first point
                                    #zd=float(z.data)
                                    zd=min(zd,Zmax) # trunkate interpolant to Zmax<=0
                                    xs.append(x[kx])
                                    ys.append(y[ky])
                                    zs.append(zd)
    if n % 10 == 0:
        print("interpolating to:"+str(yp)+":"+str(xp))
        t1 = time.time()
        time_per_iter = (t1 - t0) / (n + 1)
        time_remaining = (nn - n) * time_per_iter / 60
        print("Progress: "+str(n+1)+" of "+ str(nn))
        print(f"Time remaining: {time_remaining:.2f} minutes")
        print(f"Average time per node: {time_per_iter:.2f} seconds")
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
    
    VarObs=(np.mean(np.abs(zs))/100.)**2
#    VarObs=(np.max(np.abs(zs))/100.)**2
    #VarBG=((np.max(zs)-np.min(zs))/10.)**2
    VarBG=(np.std(zs))**2
    #stderr=np.std(zs)/10.
    if VarObs < 1.:
        VarObs  = 1.
    #VarBG=VarObs*25.
    #VarBG=1.
    
    ziID[n]=GM.InverseDistance(xs,ys,zs,xp,yp)
    ziGMU[n]=GM.GaussMarkovUnkMean(xs, ys, zs, xp, yp,LSp, VarObs,VarBG)
    ziGMN[n]=GM.GaussMarkov(xs, ys, zs, xp, yp, LSp, VarObs,VarBG,"Nearest")
    ziGMM[n]=GM.GaussMarkov(xs, ys, zs, xp, yp, LSp, VarObs,VarBG,"Mean")
    ziGM0[n]=GM.GaussMarkov(xs, ys, zs, xp, yp, LSp, VarObs,VarBG,"Zero")
    D=GM.DistanceV(xs,ys,xp,yp)
    jc=np.argmin(D)
    ziClosest[n]=zs[jc]
    

np.savetxt(floutID ,  ziID , fmt='%.6f', delimiter='\n')
np.savetxt(floutGMN,  ziGMN, fmt='%.6f', delimiter='\n')
np.savetxt(floutGMM,  ziGMM, fmt='%.6f', delimiter='\n')
np.savetxt(floutGM0,  ziGM0, fmt='%.6f', delimiter='\n')
np.savetxt(floutGMU,  ziGMU, fmt='%.6f', delimiter='\n')
np.savetxt(floutClosest,  ziClosest, fmt='%.6f', delimiter='\n')
np.savetxt(floutNpts, NumPoints, fmt='%d', delimiter='\n')
np.savetxt(floutLLS, LocalLengthScale, fmt='%.6f', delimiter='\n')
    
    
