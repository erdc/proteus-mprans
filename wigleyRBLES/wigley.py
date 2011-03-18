"""
A helper module for doing air/water flow around a moving rigid cylinder in 2D
"""
from math import *
from proteus import *

#----------------------------------------------------
# Physical properties
#----------------------------------------------------
rho_0=998.2
nu_0 =1.004e-6

rho_1=1.205
nu_1 =1.500e-5

sigma_01=0.0

g=[0.0,0.0,-9.81]

#----------------------------------------------------
# Domain - mesh - quadrature
#----------------------------------------------------
nd = 3
hull_length = 1.000
waterLevel  = 0.500 

nLevels = 1
domain = Domain.MeshTetgenDomain(fileprefix="wigley")
boundaryTags = { 'bottom': 1, 'front':2, 'right':3, 'back': 4, 'left':5, 'top':6, 'obstacle':11}

restrictFineSolutionToAllMeshes=False
parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.element
nLayersOfOverlapForParallel = 0

quad_order = 2

#----------------------------------------------------
# Boundary conditions and other flags
#----------------------------------------------------
openTop = True
openSides = True
smoothBottom = False
smoothObstacle = False
rampInitialConditions = False

movingDomain=False

checkMass=False
applyCorrection=True
applyRedistancing=True
freezeLevelSet=True

#----------------------------------------------------
# Time stepping and velocity
#----------------------------------------------------
Fr = 0.28
Um = Fr*sqrt(fabs(g[2])*hull_length)
Re = hull_length*Um*rho_0/nu_0

residence_time = hull_length/Um
dt_init=0.02
T = 5.0*residence_time

nDTout=int(ceil(T/dt_init))

print "   ================   REYNOLDS NUMBER = "+`Re`
print "   ================   FROUDE(HULL LENGTH) NUMBER = "+`Fr`
print "   ================   SPEED[M/S] = "+`Um`
print "   ================   Hull Flow time = "+`residence_time`
print "   ================   T = "+`T`
print "   ================   nDTout = "+`nDTout`

#----------------------------------------------------
# Numerical parameters
#----------------------------------------------------

useRBLES   = 0.0
useMetrics = 0.0
ns_shockCapturingFactor=0.2
ns_shockCapturingFactor=0.2

ls_shockCapturingFactor=0.2
ls_sc_uref = 1.0
ls_sc_beta = 1.5

vof_shockCapturingFactor=0.2
vof_sc_uref = 1.0
vof_sc_beta = 1.5

rd_shockCapturingFactor=0.2

epsFact_consrv_diffusion=1000.0

#----------------------------------------------------
# Interface width
#----------------------------------------------------
epsFact = 1.5

epsFact_redistance = 0.33

epsFact_density          = epsFact 
epsFact_viscosity        = epsFact 
epsFact_redistance       = epsFact 
epsFact_curvature        = epsFact 
epsFact_consrv_heaviside = epsFact 
epsFact_consrv_dirac     = epsFact 
epsFact_vof              = epsFact 



