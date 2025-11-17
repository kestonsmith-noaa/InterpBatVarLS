import os
import argparse
import time
import numpy as np
import netCDF4 as nc
import csv

import jigsawpy
import sys
from scipy.interpolate import RegularGridInterpolator
import re
#import ComputeMeshToMeshInterpWeights as mshint
import math
from geopy import distance
import itertools
#from geopy.distance import haversine
from pykrige.ok import OrdinaryKriging

import time

import FiniteElementMeshRoutines as FE
import GaussMarkov as GM

coords_1 = (34.0522, -118.2437)  # Los Angeles
coords_2 = (40.7128, -74.0060)   # New York City

#Using Haversine
#distance_haversine = distance.haversine(coords_1, coords_2).km
dist=distance.great_circle(coords_1, coords_2).km
print(f"Haversine distance: {dist:.2f} km")
maxval=100.



#[0->7] variables(dimensions): |S1 crs(), float64 lat(lat), float64 lon(lon), float32 z(lat, lon)
#[8->9]variables(dimensions): float64 x(x), float64 y(y), float32 z(y, x)
#[10]variables(dimensions): float64 lon(num_lon), float64 lat(num_lat), int16 bed_elevation(num_row, num_col)
'.S250m.nc'
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

print(nxl)
print(nyl)
print(len(xlist))

mesh="meshes/RWPS10to1km.msh"

xi, yi, ei = FE.loadWW3MeshCoords(mesh)

lsE=FE.lengthscale(xi, yi, ei)
print(lsE)

areaE=FE.ElementArea(xi, yi, ei)
print(areaE)

lsN=FE.ComputeNodeLengthScale(lsE, areaE, ei)
print(lsN)
print(np.min(lsN))
print(np.mean(lsN))
print(np.max(lsN))

nn=len(xi)
nf=len(flnms)





#nn=103
ziD=np.zeros(nn)
ziSK=np.zeros(nn)
ziOK=np.zeros(nn)
ziOKerr=np.zeros(nn)
NumPoints=np.zeros(nn)
t0 = time.time()
for n in range(nn):
    xp=xi[n]%360
    yp=yi[n]
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
    
    if n  % 200 == 0:
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
    ziOK[n]=GM.GaussMarkovUnkMean(xs, ys, zs, xp, yp,2*lsN[n], 2, stderr)
    ziSK[n]=GM.GaussMarkov(xs, ys, zs, xp, yp, 2*lsN[n], 2, stderr)

###    ziOK[n]=GM.GaussMarkovUnkMean(xs, ys, zs, xp, yp,2*lsN[n], 2, stderr)
    """
    nug=np.mean(np.abs(zs)) /10.
    sill=10.*np.mean(np.abs(zs))
    OK = OrdinaryKriging(xs, ys, zs,
        variogram_model='spherical',
        #variogram_parameters={'range': lsN[n]/111.},
    #    variogram_parameters={'psill': sill, 'range': lsN[n]/111., 'nugget': nug},
        verbose=False,
        enable_plotting=False)


    nug=np.mean(np.abs(zs)) /10.
    OK = OrdinaryKriging(xs, ys, zs,
        variogram_model='linear',
        variogram_parameters={'slope': 1000000/lsN[n], 'nugget': nug},
        verbose=False,
        enable_plotting=True)
   
#    zok, ssok = OK.execute("grid", xp,yp)
#    ziOKerr[n]=np.sqrt(ssok)
#    ziOK[n]=zok
     """
np.savetxt(mesh+"NumPoints.txt", NumPoints, fmt='%d', delimiter='\n')
np.savetxt(mesh+"DistanceBased.txt", ziD, fmt='%.6f', delimiter='\n')
np.savetxt(mesh+"SK.txt",  ziSK, fmt='%.6f', delimiter='\n')
np.savetxt(mesh+"GMOK.txt",  ziOK, fmt='%.6f', delimiter='\n')
#np.savetxt(mesh+"OK.txt",  ziOK, fmt='%.6f', delimiter='\n')
#np.savetxt(mesh+"OKstd.txt",  ziOKerr, fmt='%.6f', delimiter='\n')
    
    

"""


function lengthscale_wgs84_MEL=ComputeLengthScale_wgs84_MEL(x,y,e)
%function ComputeLengthScale(x,y,e)
% Input:
% x -lon [nn x 1]
% y -lat [nn x 1]
% e -element list[ ne x 3]
% 
% Output:
% lengthscale [ne x 1] lengthscale of element in km,
% defined by mean side length
mthd=1

x1=x(e(:,1));y1=y(e(:,1));
x2=x(e(:,2));y2=y(e(:,2));
x3=x(e(:,3));y3=y(e(:,3));

if mthd==1
    wgs84 = wgs84Ellipsoid("km");
    D3 = distance([y1,x1],[y2,x2],wgs84);
    D1 = distance([y2,x2],[y3,x3],wgs84);
    D2 = distance([y3,x3],[y1,x1],wgs84);
    lengthscale_wgs84_MEL=[D1+D2+D3]/3;% mean edge length
  %  lengthscale_wgs84_MEL=min(min(D1,D2),D3);% shortest edge
else
    R=6378.100 %radius of earth in km
    D03 = distance([y1,x1],[y2,x2],'degrees');
    D01 = distance([y2,x2],[y3,x3],'degrees');
    D02 = distance([y3,x3],[y1,x1],'degrees');
    D0=[D01+D02+D03]/3;
    lengthscale_wgs84_MEL=R*sin(D0*pi/180);
end
"""


"""
>>> z
masked_array(
  data=[[71.2074966430664, 70.46293640136719, 69.40281677246094, ..., --,
         --, --],
        [73.06659698486328, 71.84660339355469, 69.27499389648438, ...,
         --, --, --],
        [75.2576675415039, 72.89524841308594, 69.97608947753906, ..., --,
         --, --],
        ...,
        [--, --, --, ..., 37.402130126953125, 38.20088195800781,
         39.22444534301758],
        [--, --, --, ..., 36.55867385864258, 36.76195526123047,
         37.42732238769531],
        [--, --, --, ..., 35.668052673339844, 35.90686798095703,
         36.22725296020508]],
  mask=[[False, False, False, ...,  True,  True,  True],
        [False, False, False, ...,  True,  True,  True],
        [False, False, False, ...,  True,  True,  True],
        ...,
        [ True,  True,  True, ..., False, False, False],
        [ True,  True,  True, ..., False, False, False],
        [ True,  True,  True, ..., False, False, False]],
  fill_value=-99999.0,
  dtype=float32)
"""
