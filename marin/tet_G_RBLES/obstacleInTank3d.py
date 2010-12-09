"""
A helper module for doing air/water flow around a moving rigid cylinder in 2D
"""
from math import *
from proteus import *

#
#material properties
#
#water
rho_0=1000.0
nu_0=0.001/rho_0
#air
rho_1=1.0
nu_1= 0.00002/rho_1  

sigma_01=0.0
#gravity
g=[0.0,0.0,-9.81]
useStokes=False
#
#domain
#
nd = 3
height = 1.0
length = 3.22
width  = 1.0
box_height = 0.161
box_width  = 0.403
box_length = 0.161
box_xy = [2.3955,.2985]
#

useHex = False
useRBLES = 1.0
useMetrics = 1.0

if useHex:
     domain = Domain.MeshHexDomain("marinHex") 
else:
     genMesh=True
     from boxInTank3dDomain import *
     domain = boxInTank3d(L=[length,width,height],
                     box_xy=box_xy,
                     box_L=[box_length,box_width,box_height])
     domain.writePoly("boxInTank3d")
     domain.writePLY("boxInTank3d")
     domain.writeAsymptote("boxInTank3d")
     he= 3.22/64.0 
     triangleOptions="VApq1.25q12ena%f" % ((he**3)/6.0,)
     print triangleOptions

#boundaryTags = domain.boundaryTags

waterColumn_z = 0.55
waterColumn_x = 1.22

shock_x = waterColumn_x
shock_y = 10*width
shock_z = waterColumn_z

def shockSignedDistance(x):
    phi_x = x[0] - shock_x
    phi_y = x[1] - shock_y
    phi_z = x[2] - shock_z
    phi_xy = sqrt((x[0] - shock_x)**2 + (x[1] - shock_y)**2)
    phi_xz = sqrt((x[0] - shock_x)**2 + (x[2] - shock_z)**2)
    phi_yz = sqrt((x[1] - shock_y)**2 + (x[2] - shock_z)**2)
    phi_xyz = sqrt((x[0] - shock_x)**2 + (x[1] - shock_y)**2 + (x[2] - shock_z)**2)
    if x[0] > shock_x:
        if x[2] < shock_z:
            if x[1] > shock_y:
                return phi_xy
            else:
                return phi_x
        else:
            if x[1] > shock_y:
                return phi_xyz
            else:
                return phi_xz
    elif x[1] > shock_y:
        if x[2] < shock_z:
            return phi_y
        else:
            return phi_yz
    else:
        if x[2] > shock_z:
            return phi_z
        else:
            return -min(fabs(phi_x),min(fabs(phi_y),fabs(phi_z)))
#
#time interval etc.
#
T = 6.0
dt_fixed = 0.01
nDTout=int(T/dt_fixed)
dt_init=0.01
useFixedStep = True
#
#numerics
#
nLevels = 1
rdtimeIntegration='pseudoTime'#'newton'
freezeLevelSet=True

spaceOrder=1
obstacleInTank_quad_order = spaceOrder+1

useBackwardEuler=True
useBackwardEuler_ls=True
#subgrid error
lag_ns_subgridError=False
lag_ns_shockCapturing=False
lag_ls_shockCapturing=False
lag_vof_shockCapturing=False
#shock capturing diffusion
ns_shockCapturingFactor=0.2
ls_shockCapturingFactor=0.2
vof_shockCapturingFactor=0.2
rd_shockCapturingFactor=0.5
#epsilons for Heaviside/Dirac/etc smoothing
epsFact_density = 1.5
epsFact_viscosity = 1.5
epsFact_redistance = 0.33
epsFact_curvature=1.5
epsFact_consrv_heaviside=1.5
epsFact_consrv_dirac=1.5
epsFact_consrv_diffusion=10.0
epsFact_vof=1.5
#
usePETSc=True
restrictFineSolutionToAllMeshes=False
parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0
