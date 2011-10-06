"""
A helper module for doing air/water flow around a moving rigid cylinder in 2D
"""
from math import *
from proteus import *
import numpy

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
hull_length  = 5.720
hull_mass    = 532.277
hull_cg      = [2.7618104935392300,  0.0 ,0.27953462008339180  ]
hull_inertia = [[28.2823,  0.0,       20.86855 ],
                [0.0    ,  1126.799,    0.0    ],
		[20.86855,  0.0,      1136.371 ]]
			
RBR_linCons  = [1,1,0]   
RBR_angCons  = [1,0,1]  


waterLevel   = 0.241984 

nLevels = 1
domain = Domain.MeshTetgenDomain(fileprefix="mesh")
boundaryTags = { 'bottom': 1, 'front':2, 'right':3, 'back': 4, 'left':5, 'top':6, 'obstacle':7}

domainBottom = -2.5

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
dt_init=0.0025
T = 5.0*residence_time

nDTout=int(ceil(T/dt_init))

#----------------------------------------------------
# Numerical parameters
#----------------------------------------------------

useRBLES   = 0.0
useMetrics = 0.0
ns_shockCapturingFactor=0.2
ns_shockCapturingFactor=0.2

ls_shockCapturingFactor=0.1
ls_sc_uref = 1.0
ls_sc_beta = 1.5

vof_shockCapturingFactor=0.1
vof_sc_uref = 1.0
vof_sc_beta = 1.5

rd_shockCapturingFactor=0.2

epsFact_consrv_diffusion=50.0

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

#----------------------------------------------------
# Airy wave functions
#----------------------------------------------------
wave_length = 1.5     * hull_length
wave_height = 0.002   * wave_length
wave_angle  = 0.0     * pi/180.0

#----------------------------------------------------
water_depth  = waterLevel-domainBottom
wave_length  = 2.0*pi/wave_length      
wave_periode = sqrt(-g[2]*wave_length*tanh(wave_length/waterLevel))   
wave_vel_amp = wave_periode*(wave_height/sinh(wave_length*water_depth))       

xy   = lambda x:   cos(wave_angle)*x[0] + sin(wave_angle)*x[1]
kxwt = lambda x,t: wave_length*(xy(x) - Um*t) - wave_periode*t
kzh  = lambda x:   wave_length*min(x[2]-domainBottom,water_depth)

#================================================
#  Boundary conditon  lambdas
#================================================            
u_wave   = lambda x,t: wave_vel_amp * cosh(kzh(x)) * cos(kxwt(x,t)) * cos(wave_angle)  + Um  
v_wave   = lambda x,t: wave_vel_amp * cosh(kzh(x)) * cos(kxwt(x,t)) * sin(wave_angle)  
w_wave   = lambda x,t: wave_vel_amp * sinh(kzh(x)) * sin(kxwt(x,t))      
noslip   = lambda x,t: 0.0

ls_wave  = lambda x,t: -(wave_height * cos(kxwt(x,t)) + waterLevel - x[2])
vof_wave = lambda x,t: 1.0 if ls_wave(x,t) > 0.0 else 0.0

#================================================
# Print run data
#================================================
print "      Reynolds number    = "+`Re`
print "      Froude number      = "+`Fr`
print "      Hull Speed[M/S]    = "+`Um`
print "      Hull flow time[S]  = "+`residence_time`

print "      Wave length[M]     = "+`wave_length `
print "      Wave height[M]     = "+`wave_height`
print "      Wave angle[Rad]    = "+`wave_angle `
print "      Wave periode[Hz]   = "+`wave_periode`
print "      Wave velocity[M/S] = "+`wave_vel_amp `

print "      T = "+`T`
print "      nDTout = "+`nDTout`

