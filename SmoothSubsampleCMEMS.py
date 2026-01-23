import numpy as np
import netCDF4 as nc
from scipy import signal
import numpy as np

def block_mean_subsample(matrix, block_size):
    """
    Subsamples a matrix by taking the mean of each block_size x block_size block.
    """
    # Calculate new dimensions
    r, c = matrix.shape
    new_r = r // block_size
    new_c = c // block_size
    
    # Reshape into 4D array: (new_r, block_size, new_c, block_size)
    # Then take mean over axis 1 and 3 (the block dimensions)
    reshaped = matrix[:new_r*block_size, :new_c*block_size].reshape(
        new_r, block_size, new_c, block_size
    )
    NanMeanSS1=np.nanmean(reshaped,1)
    NanMeanSS12=np.nanmean(NanMeanSS1,2)#10xNX
    #return NanMeanSS1,NanMeanSS3
    return NanMeanSS12

fl="cmems_obs-sdb_glo_phy-comp_my-oa-100m-l4-s2_static.nc"
flout="cmems_obs-sdb_glo_phy-comp_my-oa-100m-l4-s2_static.PointValues500m.nc"
count=0
filval=-99999

#Get average on m *100m x m*100m grid
m=5 

data0 = nc.Dataset(fl,"r")
x=np.asarray(data0["longitude"][:])
y=np.asarray(data0["latitude"][:])
m=5
dx=x[1]-x[0]
dy=y[1]-y[0]
xp=x[0 :: m]+2*dx
yp=y[0 :: m]+2*dy
nx=len(xp)
ny=len(xp)
x0= []
y0= []
z0= []
Nblock=1000
#for ky in range(int(ny/10)):
for ky in range(ny // Nblock):
    #ky=ky+100  
    print(ky)
    print( ny//Nblock )
    pd=ky / (ny//Nblock)
    print("frac done")
    print(pd)
    ypt=yp[ky*Nblock:(ky+1)*Nblock]
    #ZT=np.asarray(data0["elevation"][ky*m*10: (ky+1)*m*10,:]) 
    ZT=np.asarray(data0["elevation"][ky*m*Nblock: (ky+1)*m*Nblock,:]) 
    ZT5=block_mean_subsample(ZT,5)
    j,k=np.where(~np.isnan(ZT5))
    xx=xp[k]
    yy=ypt[j]
    zz=ZT5[j,k]
    z0.extend(zz.tolist())
    x0.extend(xx.tolist())
    y0.extend(yy.tolist())
    print("length(z0)")
    print(len(z0))
    """
    #for kyp in range(10):
    for kyp in range(Nblock):
        nyp=ky*m+kyp
        j=np.where(~np.isnan(ZT5[kyp,:]))
        z0.extend(ZT5[kyp,j].tolist())
        x0.extend(xp[j].tolist())
        j=j[0]
        ytmp=yp[nyp]+np.zeros(len(j))
        y0.extend(ytmp.tolist())
        print("x0")
        print(x0)
        print("y0")
        print(y0)
        print("z0")
        print(z0)
    """
    np.savetxt('x0.txt', np.array(x0), fmt='%.6f', delimiter='\n')
    np.savetxt('y0.txt', np.array(y0), fmt='%.6f', delimiter='\n')
    np.savetxt('z0.txt', np.array(z0), fmt='%.6f', delimiter='\n')
    
n=len(x0)
with nc.Dataset(flout, 'w', format='NETCDF4') as ncout:
        # Create dimensions
    ncout.createDimension('points', n)  # Unlimited dimension
    
    lon_var=ncout.createVariable('lon', 'f4', ('points',))
    lon_var.units         = 'degree_east'
    lon_var.long_name     = 'longitude'
    lon_var.standard_name = 'longitude'
    lon_var.axis          = 'lon'
    lon_var[:]=x0[:]

    lat_var=ncout.createVariable('lat', 'f4', ('points',))
    lat_var.units         = 'degree_north'
    lat_var.long_name     = 'latitude'
    lat_var.standard_name = 'latitude'
    lat_var.axis          = 'lat'
    lat_var[:]=y0[:]

    z_var=ncout.createVariable('z', 'f4', ('points',),fill_value    = -99999)
#    z_var._FillValue    = -99999
    z_var.grid_mapping  = 'crs'
    z_var.long_name     = 'z'
    z_var.units         = 'meters'
    z_var.positive      = 'up'
    z_var.standard_name = 'height'
    z_var.vert_crs_name = 'EGM2008'
    z_var. vert_crs_epsg = 'EPSG:3855'
    z_var[:]=z0[:]
    
    ncout.close
