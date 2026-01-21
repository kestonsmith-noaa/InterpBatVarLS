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

import gc

"""
#too much extra resolution for our purposes
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
    "cmems_obs-sdb_glo_phy-comp_my-oa-100m-static.pts.nc",
    "crm_vol1_2023.nc.S250m.VB.nc",
    "crm_vol2_2023.nc.S250m.VB.nc",
    "crm_vol3_2023.nc.S250m.VB.nc",
    "crm_vol4_2023.nc.S250m.VB.nc",
    "crm_vol5_2023.nc.S250m.VB.nc",
    "crm_vol7_2024.nc.S250m.VB.nc",
    "crm_vol9_2023.nc.S250m.VB.nc",
    "crm_vol10_2023.nc.S250m.VB.nc",
    "crm_vol6_2023.nc.S250m.VB.nc",
    "crm_vol8_2023.nc.S250m.VB.nc",
    "crm_southak.CRMformat.nc",
    "hurl_bathy_60m_nwhi.CRMformat.nc.S250m.VB.nc",
    "PIX/ngdc_bathy_10m_wake.CRM.nc",
    "PIX/pibhmc_bathy_20m_jarvis.CRM.nc",
    "PIX/pibhmc_bathy_20m_johnston.CRM.nc",
    "PIX/pibhmc_bathy_20m_kingman.CRM.nc",
    "PIX/pibhmc_bathy_40m_baker.CRM.nc",
    "PIX/pibhmc_bathy_40m_howland.CRM.nc",
    "PIX/pibhmc_bathy_40m_palmyra.CRM.nc",
    "PIX/pibhmc_bathy_40m_rose.CRM.nc",
    "PIX/pibhmc_bathy_40m_swains.CRM.nc",
    "PIX/pibhmc_bathy_40m_vailuluu.CRM.nc",
    "PIX/pibhmc_bathy_5m_palmyra.CRM.nc",
    "PIX/sopac_bathy_50m_majuro_reef.CRM.nc",
    "PIX/ngdc_bathy_180m_mariana.CRM.nc",
    "ngdc_bathy_90m_amsamoa.crm.nc",
    "RTopo_2_0_4_GEBCO_v2023_60sec_pixel.CRMformat.nc"
    ]
nf=len(flnms)
filetype=np.zeros(nf, dtype=int)
filetype[0]=1

"""
#too many redundent data sets
flnms=[
    "crm_vol1_2023.nc.S250m.VB.nc",
    "crm_vol2_2023.nc.S250m.VB.nc",
    "crm_vol3_2023.nc.S250m.VB.nc",
    "crm_vol4_2023.nc.S250m.VB.nc",
    "crm_vol5_2023.nc.S250m.VB.nc",
    "crm_vol7_2024.nc.S250m.VB.nc",
    "crm_vol9_2023.nc.S250m.VB.nc",
    "crm_vol10_2023.nc.S250m.VB.nc", 
    "crm_vol6_2023.nc.S250m.VB.nc",
    "crm_vol8_2023.nc.S250m.VB.nc",
    "crm_southak.CRMformat.nc", # dx=700m, fillvalue= -2147483648, not used...
    "hurl_bathy_60m_nwhi.CRMformat.nc.S250m.VB.nc",
    "ngdc_bathy_90m_amsamoa.crm.nc",
    "pibhmc_bathy_60m_guam.crm.nc",
    "PIX/ngdc_bathy_10m_tutuila.CRM.nc",
    "PIX/ngdc_bathy_10m_wake.CRM.nc",
    "PIX/ngdc_bathy_180m_mariana.CRM.nc",
    "PIX/pibhmc_bathy_10m_agrihan.CRM.nc",
    "PIX/pibhmc_bathy_10m_alamagan.CRM.nc",
    "PIX/pibhmc_bathy_10m_asuncion.CRM.nc",
    "PIX/pibhmc_bathy_10m_guguan.CRM.nc",
    "PIX/pibhmc_bathy_10m_maug.CRM.nc",
    "PIX/pibhmc_bathy_10m_pagan.CRM.nc",
    "PIX/pibhmc_bathy_10m_pajaros.CRM.nc",
    "PIX/pibhmc_bathy_10m_sarigan.CRM.nc",
    "PIX/pibhmc_bathy_10m_supply.CRM.nc",
    "PIX/pibhmc_bathy_20m_jarvis.CRM.nc",
    "PIX/pibhmc_bathy_20m_johnston.CRM.nc",
    "PIX/pibhmc_bathy_20m_kingman.CRM.nc",
    "PIX/pibhmc_bathy_20m_nebank.CRM.nc",
    "PIX/pibhmc_bathy_40m_baker.CRM.nc",
    "PIX/pibhmc_bathy_40m_howland.CRM.nc",
    "PIX/pibhmc_bathy_40m_palmyra.CRM.nc",
    "PIX/pibhmc_bathy_40m_rose.CRM.nc",
    "PIX/pibhmc_bathy_40m_swains.CRM.nc",
    "PIX/pibhmc_bathy_40m_twoperbank.CRM.nc",
    "PIX/pibhmc_bathy_40m_vailuluu.CRM.nc",
    "PIX/pibhmc_bathy_5m_palmyra.CRM.nc",
    "PIX/pibhmc_bathy_5m_tinian.CRM.nc",
    "PIX/pibhmc_bathy_60m_rota.CRM.nc",
    "PIX/sopac_bathy_50m_majuro_reef.CRM.nc",
    "RTopo_2_0_4_GEBCO_v2023_60sec_pixel.CRMformat.nc"
    ]
"""
"""
flnms=[
    "crm_vol1_2023.nc.S250m.VB.nc",
    "crm_vol2_2023.nc.S250m.VB.nc",
    "crm_vol3_2023.nc.S250m.VB.nc",
    "crm_vol4_2023.nc.S250m.VB.nc",
    "crm_vol5_2023.nc.S250m.VB.nc",
    "crm_vol7_2024.nc.S250m.VB.nc",
    "crm_vol9_2023.nc.S250m.VB.nc",
    "crm_vol10_2023.nc.S250m.VB.nc", 
    "crm_vol6_2023.nc.S250m.VB.nc",
    "crm_vol8_2023.nc.S250m.VB.nc",
    "crm_southak.CRMformat.nc",# dx=700m, fillvalue= -2147483648, not used...
    "hurl_bathy_60m_nwhi.CRMformat.nc.S250m.VB.nc",
    "ngdc_bathy_90m_amsamoa.crm.nc",
    "pibhmc_bathy_60m_guam.crm.nc",
    "ngdc_bathy_180m_mariana.CRM.nc",
    "RTopo_2_0_4_GEBCO_v2023_60sec_pixel.CRMformat.nc"
    ]
"""
#    "srtm30plus_v11_bathy.CRMformat.nc"
#    "RTopo_2_0_4_GEBCO_v2023_60sec_pixel.CRMformat.nc"
#    ]
#    "crm_socal_1as_vers2.nc.S250m.nc",# mostly redundent with 6 but some extra regions, need to synthesize
Zmin=-11000. #deepest legit ocean depth value in case of mask fail 
# nescesarry for data sets "crm_vol6_2023.nc.S250m.nc" and "crm_vol8_2023.nc.S250m.nc"
Zmax=0. # maximum value to include in interpolation - zero out land values or ignore land 
Dmin=5./1000.# 5 meter minimum distance between observation points in km, prevent singularity
lambdaLL=.025 # set deg lat, lon search width for overlapping regions
#use Approximately NxTarget**2 points per data set, 
#each linear system is approx 2*(NxTarget**2)
NxTarget=20 # N=20 ~ 10 grid points in each cardinal direction used in interpolating 
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

mshnm="RWPS.OSMxGSHHS"
mesh="meshes/"+mshnm+".msh"
OutDir=mshnm+".files/"

xi, yi, ei = FE.loadWW3MeshCoords(mesh)

#lsE=FE.lengthscale(xi, yi, ei)
#areaE=FE.ElementArea(xi, yi, ei)
#lsN=FE.ComputeNodeLengthScale(lsE, areaE, ei)
lsN=FE.ComputeNodeLengthScale(xi,yi,ei)

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
floutGMUerr= OutDir+'GMU.std.'+str(N)+'.txt'
floutGMNerr= OutDir+'GMN.std.'+str(N)+'.txt'
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

#del xi yi ei lsN #clear global mesh variables
#gc.collect()

nn=len(xil)

#set search width for each data set based on expected
#number of points to interpolate
SearchWidth=np.zeros(nf)
for j in range(nf):
    x=np.array(xlist[j][:])
    dx=np.abs(x[1]-x[0])
    SearchWidth[j]=dx*float(NxTarget)/2.
    print("search width(deg Lat lon) for file: "+flnms[j]+" = ",str(SearchWidth[j]))
SearchWidth[nf-1]=.75*SearchWidth[nf-1] # smaller number for global bathy set
    
NumPoints=np.zeros(nn)
ziID=np.zeros(nn)
ziGMU=np.zeros(nn)
stdiGMU=np.zeros(nn)
ziGMN=np.zeros(nn)
stdiGMN=np.zeros(nn)
ziGMM=np.zeros(nn)
ziGM0=np.zeros(nn)
ziClosest=np.zeros(nn)
LocalLengthScale=np.zeros(nn)

t0 = time.time()
for n in range(nn):
    xp=xil[n]%360
    yp=yil[n]

#    LSp=2.*LSl[n] #probably better choice for spherical covarianve function
    LSp=1.*LSl[n] #choice for exponential covarianve function
    LSp=0.25*LSl[n] 
#    LSp=LSl[n] #Local length scale for interpolation
    LocalLengthScale[n]=LSp 
    xs=[]
    ys=[]
    zs=[]
    si=0
    for j in range(nf):
        if filetype[j] ==0: #gridded data
            if all([ xp < xmax[j],xp > xmin[j],yp < ymax[j],yp > ymin[j]]):
                x=np.array(xlist[j][:])
                y=np.array(ylist[j][:])
                fl=flnms[j]
                jx=np.array(np.where( np.abs(xp-x) < SearchWidth[j] ))
                jx=list(itertools.chain.from_iterable(jx))
                jy=np.array(np.where( np.abs(yp-y) < SearchWidth[j] ))
                jy=list(itertools.chain.from_iterable(jy))
                data = nc.Dataset(fl,"r")
                for kx in jx:
                    for ky in jy:
                            if kx and ky: # not empty 
                                z=data["z"][ky,kx]
                                if not z.mask:
                                    IncludePoint=True
                                    if len(xs)>0: #BEGIN -check for near duplicate points - causes singularity in GM
                                        d=GM.DistanceV(xs,ys,x[kx],y[ky])
                                        if np.min(d)<Dmin:
                                            print("near duplicate point at distance: "+str(np.min(d))+" for node "+str(n))
                                            IncludePoint=False #END -check for near duplicate points
                                    zd=float(z.data)
                                    if zd < Zmin:#double check mask fail and bad fill value
                                            IncludePoint=False
    #                                if zd > Zmax:#exclude land values(should be optional)
    #                                        IncludePoint=False
                                    if IncludePoint: # first point
                                        zd=min(zd,Zmax) # trunkate interpolant to Zmax<=0
                                        xs.append(x[kx])
                                        ys.append(y[ky])
                                        zs.append(zd)
        else: # global scattered data, filetype[j] == 1
            fl=flnms[j]
#            j=np.array( np.where(  np.abs(xp-x) < SearchWidth[j] and np.abs(yp-y) < SearchWidth[j]  ) )
            j=np.array( np.where( a.all( np.abs(xp-x) < SearchWidth[j], np.abs(yp-y) < SearchWidth[j] ) ) )
            j=list(itertools.chain.from_iterable(j))
            z=data["z"][j]
            zd=float(z.data)
            x=data["lon"][j]
            y=data["lat"][j]
            jp=np.where( zd <= 0.)
            x=x[jp]
            y=y[jp]
            zd=zd[jp]
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
    ObsErr= max(1. , np.mean(np.abs(zs))/100. ) # 1 percent of local mean depth
    VarObs=ObsErr**2
    #VarBG=VarObs+max(VarObs,  np.var( np.abs(zs) ) )
    #VarObs=1. 
    VarBG=10.*VarObs
    
#    print("VarObs= "+str(VarObs)+", VarBG= "+str(VarBG))
    print("std obs= " + str(np.sqrt(VarObs)) + " (m), std BG= " + str(np.sqrt(VarBG))+" (m)"  )
    
    ziID[n]=GM.InverseDistance(xs,ys,zs,xp,yp)
    ziGMU[n], stdiGMU[n]=GM.GaussMarkovUnkMean(xs, ys, zs, xp, yp,LSp, VarObs,VarBG,True)
    print("Number of obs = " +str(len(zs))  )
    print("GMU: est= "+str(ziGMU[n])+", err= "+str(stdiGMU[n]))
    ziGMN[n], stdiGMN[n]=GM.GaussMarkov(xs, ys, zs, xp, yp, LSp, VarObs,VarBG,"Nearest",True) #gives good results to eye
    print("GMN: est= "+str(ziGMN[n])+", err= "+str(stdiGMN[n]))
    ziGMM[n]=GM.GaussMarkov(xs, ys, zs, xp, yp, LSp, VarObs,VarBG,"Mean",False)
    ziGM0[n]=GM.GaussMarkov(xs, ys, zs, xp, yp, LSp, VarObs,VarBG,"Zero",False)
    D=GM.DistanceV(xs,ys,xp,yp)
    jc=np.argmin(D)
    ziClosest[n]=zs[jc]
    

np.savetxt(floutID ,  ziID , fmt='%.6f', delimiter='\n')
np.savetxt(floutGMN,  ziGMN, fmt='%.6f', delimiter='\n')
np.savetxt(floutGMNerr, stdiGMN, fmt='%.6f', delimiter='\n')
np.savetxt(floutGMM,  ziGMM, fmt='%.6f', delimiter='\n')
np.savetxt(floutGM0,  ziGM0, fmt='%.6f', delimiter='\n')
np.savetxt(floutGMU,  ziGMU, fmt='%.6f', delimiter='\n')
np.savetxt(floutGMUerr,  stdiGMU, fmt='%.6f', delimiter='\n')
np.savetxt(floutClosest,  ziClosest, fmt='%.6f', delimiter='\n')
np.savetxt(floutNpts, NumPoints, fmt='%i', delimiter='\n')
np.savetxt(floutLLS, LocalLengthScale, fmt='%.6f', delimiter='\n')
    
    
