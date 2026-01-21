
import FiniteElementMeshRoutines as FE
import numpy as np
import sys

mshnm="RWPS.OSMxGSHHS"
mesh="meshes/"+mshnm+".msh"

flin=sys.argv[1]
flout=flin+".WW3.msh"

x, y, z0, e, bnd = FE.loadWW3Mesh(mesh)

zmin=1.
nn=len(x)
z=np.zeros(nn)
f=open(flin, 'r')
for k in range(nn):
    line=f.readline()
    #print(str(k)+" "+line)
    z[k]=float(line.strip())
    z[k]=-z[k]
    z[k]=max(z[k],zmin)
f.close
FE.WriteWW3Mesh(flout,x,y,z,e,bnd)

