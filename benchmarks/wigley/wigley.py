"""
A helper module for doing air/water flow around a moving rigid cylinder in 2D
"""
from math import *
from proteus import *
import numpy
from proteus.ctransportCoefficients import smoothedHeaviside

#----------------------------------------------------
# Physical properties
#----------------------------------------------------
rho_0=998.2
nu_0 =1.004e-6

rho_1=rho_0 #1.205
nu_1 =1.500e-5

sigma_01=0.0

g=[0.0,0.0,-9.81]

#----------------------------------------------------
# Domain - mesh - quadrature
#----------------------------------------------------
nd = 3
hull_length  = 1.0
hull_mass    = 532.277
hull_cg      = [2.7618104935392300,  0.0 ,0.27953462008339180  ]
hull_inertia = [[28.2823,  0.0,       20.86855 ],
                [0.0    ,  1126.799,    0.0    ],
		[20.86855,  0.0,      1136.371 ]]
			
RBR_linCons  = [1,1,0]   
RBR_angCons  = [1,0,1]  


waterLevel   = 0.5

nLevels = 1
#domain = Domain.MeshTetgenDomain(fileprefix="mesh")
#boundaryTags = { 'bottom': 1, 'front':2, 'right':3, 'back': 4, 'left':5, 'top':6, 'obstacle':11}
L=(2.75,2.0, 0.8)
Refinement = 4
if False:
  nnx=16
  nny=8
  nnz=16
  L=(2.75,2.0, 0.8)
  domain = Domain.RectangularDomain()#L=[2.75,2.0, 0.8])
  nx = 40
  he = 2.75/(nx-1)
  triangleOptions="VApq2q10ena%21.16e" % ((he**3)/6.0,)
  print triangleOptions
  domain.writePoly("mesh")
  boundaryTags = domain.boundaryLegend
  boundaryTags['obstacle'] = 11 

else:

  he = L[0]/float(4*Refinement-1)
  boundaries=['left','right','bottom','top','front','back','obstacle']
  boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
  vertices=[[0.0,0.0,0.0],#0
	    [L[0],0.0,0.0],#1
	    [L[0],L[1],0.0],#2
	    [0.0,L[1],0.0],#3
	    [0.0,0.0,L[2]],#4
	    [L[0],0.0,L[2]],#5
	    [L[0],L[1],L[2]],#6
	    [0.0,L[1],L[2]]]#7
  vertexFlags=[boundaryTags['left'],
	     boundaryTags['right'],
	     boundaryTags['right'],
	     boundaryTags['left'],
	     boundaryTags['left'],
	     boundaryTags['right'],
	     boundaryTags['right'],
	     boundaryTags['left']]
  facets=[[[0,1,2,3]],
  	  [[0,1,5,4]],
	  [[1,2,6,5]],
	  [[2,3,7,6]],
	  [[3,0,4,7]],
	  [[4,5,6,7]]]
  facetFlags=[boundaryTags['bottom'],
	      boundaryTags['front'],
	      boundaryTags['right'],
	      boundaryTags['back'],
	      boundaryTags['left'],
	      boundaryTags['top']]
  regions=[[0.5*L[0],0.5*L[1],0.5*L[2]]]
  regionFlags=[1.0]
  domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
					       vertexFlags=vertexFlags,
					       facets=facets,
					       facetFlags=facetFlags,
					       regions=regions,
					       regionFlags=regionFlags)
  #go ahead and add a boundary tags member 
  domain.boundaryTags = boundaryTags
  domain.writePoly("mesh")
  domain.writePLY("mesh")
  domain.writeAsymptote("mesh")
  triangleOptions="VApq2q10ena%21.16e" % ((he**3)/6.0,)






domainBottom = -2.5

restrictFineSolutionToAllMeshes=False
parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.element
nLayersOfOverlapForParallel = 0
use_petsc4py=True#Original PETSc solvers do not appear to be working: 01/10/13 

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
dt_init=0.025
T = 5.0*residence_time

nDTout=int(ceil(T/dt_init))

#----------------------------------------------------
# Numerical parameters
#----------------------------------------------------

useRBLES   = 0.0
useMetrics = 0.0
ns_shockCapturingFactor=0.2
ns_shockCapturingFactor=0.2

ls_shockCapturingFactor=0.05
ls_sc_uref = 1.0
ls_sc_beta = 1.0

vof_shockCapturingFactor=0.05
vof_sc_uref = 1.0
vof_sc_beta = 1.0

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
wave_length = 1.5e8   * hull_length
wave_height = 0.00    * wave_length
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

u_wave = lambda x,t : Um
v_wave = lambda x,t : 0.0
w_wave = v_wave

noslip   = lambda x,t: 0.0
noflow   = lambda x,t: 0.0

hs_pres  = lambda x,t: -g[2]*max(rho_0*(waterLevel-x[2]),rho_1*(waterLevel-x[2]))

ls_wave  = lambda x,t: -(wave_height * cos(kxwt(x,t)) + waterLevel - x[2])
#vof_wave = lambda x,t: 1.0 if ls_wave(x,t) > 0.0 else 0.0

vof_wave = lambda x,t: smoothedHeaviside(epsFact_density*0.3,ls_wave(x,t) )

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

