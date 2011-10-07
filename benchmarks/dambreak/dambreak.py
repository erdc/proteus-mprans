from math import *
import proteus.MeshTools
from proteus import Domain
from proteus.default_n import *   
   
#  Discretization -- input options  
Refinement=5
spaceOrder=1
useHex=False
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
L = (0.584,0.146,0.350)

nLevels = 1
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
nLayersOfOverlapForParallel = 0

if useHex:   
    nnx=4*Refinement
    nny=1*Refinement
    nnz=2*Refinement
    hex=True    
    domain = Domain.RectangularDomain(L)

else:
    he = L[0]/float(4*Refinement-1)
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


# Time stepping
T=0.40
dt_fixed = 0.04/Refinement
nDTout = int(round(T/dt_fixed))

# Numerical parameters
ns_shockCapturingFactor  = 0.2
ls_shockCapturingFactor  = 0.2
ls_sc_uref  = 1.0
ls_sc_beta  = 1.0
vof_shockCapturingFactor = 0.2
vof_sc_uref = 1.0
vof_sc_beta = 1.0
rd_shockCapturingFactor  = 0.2

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

# Initial condition
waterLine_x = 0.146
waterLine_z = 0.292

def signedDistance(x):
    phi_x = x[0]-waterLine_x
    phi_z = x[2]-waterLine_z 
    if phi_x < 0.0:
        if phi_z < 0.0:
            return max(phi_x,phi_z)
        else:
            return phi_z
    else:
        if phi_z < 0.0:
            return phi_x
        else:
            return sqrt(phi_x**2 + phi_z**2)



