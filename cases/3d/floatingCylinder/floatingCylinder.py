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
#rho_1=rho_0
#nu_1=nu_0
sigma_01=0.0#72.8e-3
#gravity
g=[0.0,0.0,-9.8]
useStokes=False
#
#domain
#
nd = 3
inflow_height=0.41
he = (0.25)*inflow_height
#he = (4.0*0.025)*inflow_height
#he = (3.0*0.025)*inflow_height
#he = 0.025*inflow_height
inflow_width=he#4.0*he
bottom_length=2.2
cylinder_radius=inflow_height/5.0
cylinder_center = (bottom_length/2.0,inflow_width/2.0,inflow_height/2.0)#(1.0*bottom_length/4.0,inflow_width/2.0,inflow_height/2.0)
from cylinder3dDomain import *
domain = cylinder3d(fileprefix="cylinder",
                    cross_section=circular_cross_section,
                    height=inflow_height,
                    width=inflow_width,
                    length=bottom_length,
                    radius=cylinder_radius,
                    center = cylinder_center,
#                    n_points_on_obstacle=2*5-2)
                    n_points_on_obstacle=2*21-2)
domain.writePoly("cylinder")
domain.writePLY("cylinder")
boundaryTags = domain.boundaryTags
print boundaryTags
waterLevel = 0.5*inflow_height
useShock=False#True
movingDomain=False#True
altBC=True

#
#residence time based on mean velocity
#
RE = 100000.0#100000.0
print "RE",RE
Um = nu_0*RE/(2.0*cylinder_radius)
#Um = 2.0
RE = 2.0*cylinder_radius*Um/nu_0
print "RE",RE
Profiling.logEvent("REYNOLDS NUMBER = "+`RE`)
residence_time = bottom_length/Um
#Um=0.0
#
#time interval etc.
#
dt_init=he*0.001#0.0001*residence_time
T = 1.0*residence_time
print "T = ",T
nDTout=400#1000
runCFL = 0.33
#
#numerics
#
nLevels = 1
triangleOptions="VApq1.15q15ena%21.16e" % ((he**3)/6.0,)
#triangleOptions="pAq30Dena%f" % (0.5*he**2,)
print triangleOptions
applyCorrection=True
applyRedistancing=True
#rdtimeIntegration='osher'
rdtimeIntegration='newton'
freezeLevelSet=True
quad_order = 3
useBackwardEuler=True
useBackwardEuler_ls=True
#subgrid error
lag_ns_subgridError=True
lag_ns_shockCapturing=True
lag_ls_shockCapturing=True
#shock capturing diffusion
ns_shockCapturingFactor=0.9
ls_shockCapturingFactor=0.9
vof_shockCapturingFactor=0.9
rd_shockCapturingFactor=0.9
#epsilons for Heaviside/Dirac/etc smoothing
hFactor=1.0
noPressureStabilization=False
epsFact_density = 1.5
epsFact_viscosity = 1.5
epsFact_redistance = 0.33
epsFact_curvature=1.5
epsFact_consrv_heaviside=1.5
epsFact_consrv_dirac=1.5
epsFact_consrv_diffusion=10.0
useConstantH=True
#
usePETSc=False
spaceOrder=1
restrictFineSolutionToAllMeshes=False
parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 1
