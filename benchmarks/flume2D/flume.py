from math import *
import proteus.MeshTools
from proteus import Domain
from proteus.default_n import *   
   
#  Discretization -- input options  
genMesh=True
useOldPETSc=False
useSuperlu=True
spaceOrder = 1
useHex     = False
useRBLES   = 0.0
useMetrics = 1.0
applyCorrection=True
useVF = 1.0
useOnlyVF = False
redist_Newton = False#True
useRANS = 1 # 0 -- None
            # 1 -- K-Epsilon
            # 2 -- K-Omega
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
nd = 2
if spaceOrder == 1:
    hFactor=1.0
    if useHex:
	 basis=C0_AffineLinearOnCubeWithNodalBasis
         elementQuadrature = CubeGaussQuadrature(nd,2)
         elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,2)     	 
    else:
    	 basis=C0_AffineLinearOnSimplexWithNodalBasis
         elementQuadrature = SimplexGaussQuadrature(nd,3)
         elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3) 	    
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
#L = (8.0,1.0)
L = (2.0,1.0)
he = L[1]/10
he*=0.5
he*=0.5
#print he
useObstacle=True#False
nLevels = 1
#parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0
structured = False#True
if useHex:   
    nnx=4*Refinement+1
    nny=2*Refinement+1
    hex=True    
    domain = Domain.RectangularDomain(L)
else:
    boundaries=['left','right','bottom','top','front','back','obstacle']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    obstacle_radius=0.1*L[1]
    if structured:
        nnx = int(round(L[0]/he) + 1)
        nny = int(round(L[1]/he) + 1)
    else:
        if useObstacle:
            vertices=[[0.0,0.0],#0
                      [L[0],0.0],#1
                      [L[0],L[1]],#2
                      [0.0,L[1]],#3
                      [0.5*L[0]-obstacle_radius,0.5*L[1]-obstacle_radius],#4
                      [0.5*L[0]+obstacle_radius,0.5*L[1]-obstacle_radius],#5
                      [0.5*L[0]+obstacle_radius,0.5*L[1]+obstacle_radius],#6
                      [0.5*L[0]-obstacle_radius,0.5*L[1]+obstacle_radius]]#7
            vertexFlags=[boundaryTags['bottom'],
                         boundaryTags['bottom'],
                         boundaryTags['top'],
                         boundaryTags['top'],
                         boundaryTags['obstacle'],
                         boundaryTags['obstacle'],
                         boundaryTags['obstacle'],
                         boundaryTags['obstacle']]
            segments=[[0,1],
                      [1,2],
                      [2,3],
                      [3,0],
                      [4,5],
                      [5,6],
                      [6,7],
                      [7,4]]
            segmentFlags=[boundaryTags['bottom'],
                          boundaryTags['right'],
                          boundaryTags['top'],
                          boundaryTags['left'],
                          boundaryTags['obstacle'],
                          boundaryTags['obstacle'],
                          boundaryTags['obstacle'],
                          boundaryTags['obstacle']]
            regions=[[0.5*L[0],0.5*L[1]]]
            regionFlags=[1]
            holes = [[0.5*L[0],0.5*L[1]]]
            domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                                          vertexFlags=vertexFlags,
                                                          segments=segments,
                                                          segmentFlags=segmentFlags,
                                                          regions=regions,
                                                          regionFlags=regionFlags,
                                                          holes=holes)
        else:
            vertices=[[0.0,0.0],#0
                      [L[0],0.0],#1
                      [L[0],L[1]],#2
                      [0.0,L[1]]]#3
            vertexFlags=[boundaryTags['bottom'],
                         boundaryTags['bottom'],
                         boundaryTags['top'],
                         boundaryTags['top']]
            segments=[[0,1],
                      [1,2],
                      [2,3],
                      [3,0]]
            segmentFlags=[boundaryTags['bottom'],
                          boundaryTags['right'],
                          boundaryTags['top'],
                          boundaryTags['left']]
            regions=[[0.5*L[0],0.5*L[1]]]
            regionFlags=[1]
            domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                                          vertexFlags=vertexFlags,
                                                          segments=segments,
                                                          segmentFlags=segmentFlags,
                                                          regions=regions,
                                                          regionFlags=regionFlags)       
        #go ahead and add a boundary tags member 
        domain.boundaryTags = boundaryTags
        domain.writePoly("mesh")
        domain.writePLY("mesh")
        domain.writeAsymptote("mesh")
        triangleOptions="VApq30Dena%8.8f" % ((he**2)/2.0,)

# Numerical parameters
ns_forceStrongDirichlet = False
if useMetrics:
    ns_shockCapturingFactor  = 0.1
    ns_lag_shockCapturing = True#False
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.1
    ls_lag_shockCapturing = True#False
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.0
    vof_shockCapturingFactor = 0.1
    vof_lag_shockCapturing = True#False
    vof_sc_uref = 1.0
    vof_sc_beta = 1.0
    rd_shockCapturingFactor  = 0.9
    rd_lag_shockCapturing = False
    epsFact_density    = 1.5
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 10.0
    redist_Newton = False
    kappa_shockCapturingFactor = 0.1
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.1
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.0
else:
    ns_shockCapturingFactor  = 0.9
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.9
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.0
    vof_shockCapturingFactor = 0.9
    vof_lag_shockCapturing = True
    vof_sc_uref  = 1.0
    vof_sc_beta  = 1.0
    rd_shockCapturingFactor  = 0.9
    rd_lag_shockCapturing = False
    epsFact_density    = 1.5
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 10.0
    redist_Newton = False
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref  = 1.0
    kappa_sc_beta  = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref  = 1.0
    dissipation_sc_beta  = 1.0

ns_nl_atol_res = max(1.0e-8,0.1*he**2/2.0)
vof_nl_atol_res = max(1.0e-8,0.1*he**2/2.0)
ls_nl_atol_res = max(1.0e-8,0.1*he**2/2.0)
rd_nl_atol_res = max(1.0e-8,0.1*he)
mcorr_nl_atol_res = max(1.0e-8,0.01*he**2/2.0)
kappa_nl_atol_res = max(1.0e-8,0.01*he**2/2.0)
dissipation_nl_atol_res = max(1.0e-8,0.01*he**2/2.0)

#turbulence
ns_closure=0 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
if useRANS == 1:
    ns_closure = 3
elif useRANS == 2:
    ns_closure == 4
# Water
rho_0 = 998.2
nu_0  = 1.004e-6

# Air
rho_1 = 1.205
nu_1  = 1.500e-5 

# Surface tension
sigma_01 = 0.0

# Gravity
g = [0.0,-9.8]


from proteus.ctransportCoefficients import smoothedHeaviside_integral
from proteus.ctransportCoefficients import smoothedHeaviside

inflowHeight = L[1]/2.0

Fr = 2.0#1.25
#Fr = 0.25
Um = Fr*sqrt(fabs(g[1])*2*obstacle_radius)
Re = 2*obstacle_radius*Um/nu_0
Frd = Um/sqrt(fabs(g[1])*inflowHeight)
Re_dw = Um*inflowHeight/nu_0
Re_da = Um*inflowHeight/nu_1
#print Fr,Re,Um
residence_time = L[0]/Um
inflowVelocity = (Um,0.0)
outflowVelocity = (Um,0.0)
#for RANS
kInflow = 0.003*Um

# Time stepping
T=1.0*residence_time#5.0*residence_time
dt_fixed = residence_time/40.0
dt_init = min(0.1*dt_fixed,0.001)
runCFL=0.33#33
nDTout = int(round(T/dt_fixed))
def signedDistance(x):
    return x[1]-inflowHeight

def inflowPhi(x,t):
    return signedDistance(x)

def inflowVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,x[1] - inflowHeight)

def inflowPressure(x,t):
    p_L = L[1]*rho_1*g[1]
    phi_L = L[1] - inflowHeight
    phi = x[1] - inflowHeight
    return p_L -g[1]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_density*he,phi_L)
                                                          -smoothedHeaviside_integral(epsFact_density*he,phi)))

outflowHeight = L[1]/2.0

def outflowPhi(x,t):
    return signedDistance(x)

def outflowVF(x,t):
    return smoothedHeaviside(epsFact_density*he,x[1] - outflowHeight)

def outflowPressure(x,t):
    p_L = L[1]*rho_1*g[1]
    phi_L = L[1] - outflowHeight
    phi = x[1] - outflowHeight
    return p_L -g[1]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_density*he,phi_L)
                                                          -smoothedHeaviside_integral(epsFact_density*he,phi)))

#from pylab import *
#z = np.linspace(0.0,L[1],512)
#p = np.array([outflowPressure([L[0],zi],0.0) for zi in z])
#p2 = np.array([inflowPressure([L[0],zi],0.0) for zi in z])
#plot(z,p)
#plot(z,p2)
#show()
logEvent("""
Reynolds number    = %16.21e
Froude number(obj) = %16.21e
Froude number(dep) = %16.21e
Reynold number_dw  = %16.21e
Reynold number_da  = %16.21e
Hull Speed[M/S]    = %16.21e
Hull flow time[S]  = %16.21e

T                  = %16.21e
nDTout             = %i

he                 = %16.21e
""" % (Re,           
       Fr,       
       Frd,
       Re_dw,
       Re_da,
       Um,           
       residence_time,
       T,
       nDTout,
       he))
