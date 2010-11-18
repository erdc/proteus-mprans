"""
A helper module for doing air/water flow around a moving rigid cylinder in 2D
"""
from math import *
from proteus import *

#
#material properties
#
#water
rho_0=998.2
nu_0=1.004e-6
#air
rho_1=1.205
nu_1= 1.500e-5
#rho_0=rho_1
#nu_0=nu_1
sigma_01=0.0#72.8e-3
#gravity
g=[0.0,0.0,-9.8]
useStokes=False
#
#domain
#
nd = 3
inflow_height=0.41
he = (4.0*0.025)*inflow_height
#he = (2.0*0.025)*inflow_height
#he = 0.025*inflow_height
inflow_width=he
bottom_length=0.82#2.2
cylinder_radius=inflow_height/10.0
cylinder_center = (1.0*bottom_length/4.0,inflow_width/2.0,inflow_height/2.0)
from cylinder3dDomain import *
domain = cylinder3d(fileprefix="cylinder",
                    cross_section=circular_cross_section,
                    height=inflow_height,
                    width=inflow_width,
                    length=bottom_length,
                    radius=cylinder_radius,
                    center = cylinder_center,
                    n_points_on_obstacle=2*11-2)
domain.writePoly("cylinder")
domain.writePLY("cylinder")
boundaryTags = domain.boundaryTags
waterLevel = 0.5*inflow_height
useShock=False#True
movingDomain=False#True
#
#residence time based on mean velocity
#
#RE = 100.0
Um = 0.5#nu_0*RE/(2.0*cylinder_radius)
RE = 2.0*cylinder_radius*Um/nu_0
Profiling.logEvent("REYNOLDS NUMBER = "+`RE`)
residence_time = bottom_length/Um
#Um=0.0
#
#time interval etc.
#
dt_init=0.01#0.001*residence_time
T = 1.0#10.0*residence_time
nDTout=100
runCFL = 0.33
#
#numerics
#
nLevels = 1
triangleOptions="VApq1.25q12ena%21.16e" % ((he**3)/6.0,)
#triangleOptions="pAq30Dena%f" % (0.5*he**2,)
print triangleOptions
applyCorrection=True
applyRedistancing=True
rdtimeIntegration='osher'
freezeLevelSet=True
quad_order = 3
useBackwardEuler=True
useBackwardEuler_ls=True
#subgrid error
lag_ns_subgridError=True
lag_ns_shockCapturing=True
lag_ls_shockCapturing=True
#shock capturing diffusion
ns_shockCapturingFactor=0.33
ls_shockCapturingFactor=0.33
vof_shockCapturingFactor=0.33
rd_shockCapturingFactor=0.33
#epsilons for Heaviside/Dirac/etc smoothing
hFactor=1.0
noPressureStabilization=False
epsFact_density = 1.5
epsFact_viscosity = 1.5
epsFact_redistance = 0.75
epsFact_curvature=1.5
epsFact_consrv_heaviside=1.5
epsFact_consrv_dirac=1.5
epsFact_consrv_diffusion=10.0
epsFact_vof=1.5
#
usePETSc=False#True
spaceOrder=1
restrictFineSolutionToAllMeshes=False
parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 1
