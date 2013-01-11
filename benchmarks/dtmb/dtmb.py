"""
A helper module for doing air/water flow around a moving rigid cylinder in 2D
"""
from math import *
from proteus import *
import numpy
#----wavetank info-----
from math import *
import proteus.MeshTools
import numpy as np
import waveModules_Matt as wm
from proteus import Domain
from proteus.Profiling import logEvent
from proteus.default_n import *   
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral

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
L = ()
restrictFineSolutionToAllMeshes=False
parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.node
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
dt_init=0.000025
T = 5.0*residence_time
nDTout=500
dt_out =  (T-dt_init)/nDTout
runCFL = 0.33
he = 5.0

#----------------------------------------------------
# Numerical parameters
#----------------------------------------------------

useRBLES   = 0.0
useMetrics = 0.0
ns_shockCapturingFactor=0.9

ls_shockCapturingFactor=0.9
ls_sc_uref = 1.0
ls_sc_beta = 1.5

vof_shockCapturingFactor=0.9
vof_sc_uref = 1.0
vof_sc_beta = 1.5

rd_shockCapturingFactor=0.9

epsFact_consrv_diffusion=10.0


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
logEvent("""
Reynolds number    = %16.21e
Froude number      = %16.21e
Hull Speed[M/S]    = %16.21e
Hull flow time[S]  = %16.21e

Wave length[M]     = %16.21e
Wave height[M]     = %16.21e
Wave angle[Rad]    = %16.21e
Wave periode[Hz]   = %16.21e
Wave velocity[M/S] = %16.21e
T                  = %16.21e
nDTout             = %i
""" % (Re,           
       Fr,           
       Um,           
       residence_time,
       wave_length,
       wave_height,  
       wave_angle,  
       wave_periode, 
       wave_vel_amp,
       T,
       nDTout))

#  Discretization -- input options  
Refinement = 1#4#15
genMesh=False
useOldPETSc=False
useSuperlu = False # set to False if running in parallel with petsc.options
spaceOrder = 1
useHex     = False
useRBLES   = 0.0
useMetrics = 0.0

# Input checks
if spaceOrder not in [1,2]:
    print "INVALID: spaceOrder" + spaceOrder
    sys.exit()    
    
if useRBLES not in [0.0, 1.0]:
    print "INVALID: useRBLES" + useRBLES 
    sys.exit()

if useMetrics not in [0.0, 1.0]:
    print "INVALID: useMetrics"
    sys.exit()
    
#  Discretization   
nd = 3
if spaceOrder == 1:
    hFactor=1.0
    if useHex:
	 basis=C0_AffineLinearOnCubeWithNodalBasis
         elementQuadrature = CubeGaussQuadrature(nd,2)
         elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,2)     	 
    else:
    	 basis=C0_AffineLinearOnSimplexWithNodalBasis
         elementQuadrature = SimplexGaussQuadrature(nd,2)
         elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,2) 	    
elif spaceOrder == 2:
    hFactor=0.5
    if useHex:    
	basis=C0_AffineLagrangeOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd,4)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,4)    
    else:    
	basis=C0_AffineQuadraticOnSimplexWithNodalBasis	
        elementQuadrature = SimplexGaussQuadrature(nd,4)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)
    

#wave/current properties
windspeed_u = Um
windspeed_v = 0.0
windspeed_w = 0.0

outflowHeight = waterLevel
outflowVelocity = (0.0,0.0,0.0)#not used for now

inflowHeightMean = waterLevel#0.5*L[2]
inflowVelocityMean = (0.0,0.0,0.0)

rightEndClosed=False

def waveHeight(x,t):
    return waterLevel

def waveVelocity_u(x,t):
    return Um

def waveVelocity_v(x,t):
    return 0.0

def waveVelocity_w(x,t):
    return 0.0

def wavePhi(x,t):
    return x[2] - waveHeight(x,t)

def wavePhi_init(x,t):
    return wavePhi(x,t)

def waveVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t))

def waveVF_init(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi_init(x,t))

def twpflowVelocity_u(x,t):
    waterspeed = waveVelocity_u(x,t)
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
    return H*windspeed_u + (1.0-H)*waterspeed

def twpflowVelocity_v(x,t):
    waterspeed = waveVelocity_v(x,t)
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
    return H*windspeed_v+(1.0-H)*waterspeed

def twpflowVelocity_w(x,t):
    waterspeed = waveVelocity_w(x,t)
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
    return H*windspeed_w+(1.0-H)*waterspeed

def twpflowVelocity_u_init(x,t):
    #return 0.0 # for flat initial mean water surface
    return twpflowVelocity_u(x,t)

def twpflowVelocity_v_init(x,t):
    #return 0.0 # for flat initial surface
    return twpflowVelocity_v(x,t)

def twpflowVelocity_w_init(x,t):
    #return 0.0 # for flat initial surface
    return twpflowVelocity_w(x,t)

def twpflowFlux(x,t):
    return -twpflowVelocity_u(x,t)

def outflowVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,x[2] - outflowHeight)

Lz = 2.5 - (-3.25)

def twpflowPressure(x,t):
    p_L = Lz*rho_1*g[2]
    phi_L = wavePhi((x[0],x[1],Lz),t) 
    phi = wavePhi(x,t)
    return p_L - g[2]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi_L)
                                                          -smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi)))
def twpflowPressure_init(x,t):
    p_L = Lz*rho_1*g[2]
    phi_L = Lz - inflowHeightMean
    phi = x[2] - inflowHeightMean
    return p_L -g[2]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi_L)
                                                          -smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi)))

def outflowPressure(x,t):
    p_L = Lz*rho_1*g[2]
    phi_L = Lz - outflowHeight
    phi = x[2] - outflowHeight
    return p_L -g[2]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi_L)
                                                          -smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi)))


