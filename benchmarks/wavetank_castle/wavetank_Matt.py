from math import *
import proteus.MeshTools
import waveModules as wm
from proteus import Domain
from proteus.default_n import *   
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
   
#  Discretization -- input options  
Refinement = 3 #15
genMesh=True
useOldPETSc=False
useSuperlu = True
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
    
# Domain and mesh
L = (10.0,
     10.0,
     1.0)

quasi2D = True

nLevels = 1
#parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0

if useHex:   
    nnx=4*Refinement
    nny=1*Refinement
    nnz=2*Refinement
    hex=True    
    domain = Domain.RectangularDomain(L)

else:
    he = L[2]/float(4*Refinement-1)
    if quasi2D:
        L = (L[0],he,L[2])
    boundaries=['left','right','bottom','top','front','back']
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


# Numerical parameters
ns_shockCapturingFactor  = 0.9
ls_shockCapturingFactor  = 0.9
ls_sc_uref  = 1.0
ls_sc_beta  = 1.0
vof_shockCapturingFactor = 0.9
vof_sc_uref = 1.0
vof_sc_beta = 1.0
rd_shockCapturingFactor  = 0.9

epsFact_density    = 1.5
epsFact_viscosity  = 1.5
epsFact_redistance = 0.33
epsFact_curvature  = 1.5
epsFact_consrv_heaviside = 1.5
epsFact_consrv_dirac     = 1.5
epsFact_consrv_diffusion = 10.0
epsFact_vof = 1.5

# Water
rho_0 = 998.2
nu_0  = 1.004e-6

# Air
rho_1 = 1.205
nu_1  = 1.500e-5 

# Surface tension
sigma_01 = 0.0

# Gravity
g = [0.0,0.0,-9.8]

#wave/current properties
windspeed_u = 0.0
windspeed_v = 0.0
windspeed_w = 0.0

outflowHeight = 0.5*L[2]
outflowVelocity = (0.0,0.0,0.0)#not used for now

inflowHeightMean = 0.5*L[2]
inflowVelocityMean = (0.2,0.0,0.0)

waveLength = inflowHeightMean*5 #
period = waveLength/sqrt((-g[2])*inflowHeightMean) #meters
omega = 2.0*pi/period
k=2.0*pi/waveLength
amplitude = 0.1*inflowHeightMean

def waveHeight(x,t):
    z = x[2] - inflowHeightMean
#    return inflowHeightMean + amplitude*sin(omega*t-k*x[0])
#    return inflowHeightMean + wm.initialWaveHeight(amplitude, omega, k, t, x[0])
    return inflowHeightMean + wm.linear1DHeight(amplutude, omega, k, t, x[0], depth=L[1], z)

def waveVelocity_u(x,t):
#    return inflowVelocityMean[0] + omega*amplitude*sin(omega*t - k*x[0])/(k*inflowHeightMean)
    z = x[2] - inflowHeightMean
    return inflowVelocityMean[0] + wm.linear1DVelocity(amplutude, omega, k, t, x[0], depth=L[1],z)

def waveVelocity_w(x,t):
    z = x[2] - inflowHeightMean
    return inflowVelocityMean[2] + (z + inflowHeightMean)*omega*amplitude*cos(omega*t-k*x[0])/inflowHeightMean

####
def wavePhi(x,t):
    return x[2] - waveHeight(x,t)

def waveVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t))


def twpflowVelocity_u(x,t):
    waterspeed = waveVelocity_u(x,t)
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
    return H*windspeed_u + (1.0-H)*waterspeed

def twpflowVelocity_v(x,t):
    waterspeed = 0.0
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
    return H*windspeed_v+(1.0-H)*waterspeed

def twpflowVelocity_w(x,t):
    waterspeed = waveVelocity_w(x,t)
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
    return H*windspeed_w+(1.0-H)*waterspeed

def twpflowFlux(x,t):
    return -twpflowVelocity_u(x,t)

def outflowVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,x[2] - outflowHeight)

def twpflowPressure(x,t):
    p_L = L[2]*rho_1*g[2]
    phi_L = wavePhi((x[0],x[1],L[2]),t) 
    phi = wavePhi(x,t)
    return p_L -g[2]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi_L)
                                                          -smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi)))

def outflowPressure(x,t):
    return twpflowPressure((0.0,0.0,x[2]),0.0)

# Time 
T=period*100
print "T",T
dt_fixed = period/5.0 
dt_init = 1.0e-4
nDTout = int(T/dt_fixed)
tnList = [i*dt_fixed for i in range(0,nDTout+1)] 
print tnList
