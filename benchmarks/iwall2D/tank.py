from math import *
import proteus.MeshTools
from proteus import Domain
from proteus.default_n import *   
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
import numpy as np

transect = np.loadtxt('transectVeryShort.txt',skiprows=0,delimiter="\t")
benchStart=2
benchEnd=7
# transect = np.loadtxt('transectShort.txt',skiprows=0,delimiter="\t")
# benchStart=15
# benchEnd=20
##transect = np.loadtxt('transect.txt',skiprows=0,delimiter="\t")
#transect = np.loadtxt('transect_001_bathy_meters.csv',skiprows=1,delimiter="\t")
#wse = np.loadtxt('wse_by_wavemaker_storm36_wave3.csv',skiprows=2,delimiter="\t")

#wave generator
windVelocity = (0.0,0.0)
inflowHeightMean = 2.64
inflowVelocityMean = (0.0,0.0)
period = 2.0
omega = 2.0*math.pi/period
waveheight = 0.4
amplitude = waveheight/ 2.0
wavelength = 6.2
k = 2.0*math.pi/wavelength

#  Discretization -- input options  
genMesh=True
useOldPETSc=False
useSuperlu=True
timeDiscretization='vbdf'#'vbdf'#'be','flcbdf'
spaceOrder = 1
useHex     = False
useRBLES   = 0.0
useMetrics = 1.0
applyCorrection=True
useVF = 1.0
useOnlyVF = False
useRANS = 0 # 0 -- None
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
#L = (40.0,0.7)
#for debugging, make the tank short
#L = (10.0,0.7)
#he = L[1]/10.0 #try this first
#he*=0.5
#he*=0.5
#he*=0.5
nLevels = 1
#parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0
structured=False
if useHex:   
    nnx=ceil(L[0]/he)+1
    nny=ceil(L[1]/he)+1
    hex=True    
    domain = Domain.RectangularDomain(L)
else:
    maxZ = max(transect[:,1])
    minZ = min(transect[:,1])
    maxX = max(transect[:,0])
    minX = min(transect[:,0])
    Lx = maxX - minX
    Lz = 1.3*(maxZ-minZ)
    L = (Lx,Lz,1.0)
    coarse_he = wavelength*0.1
    #coarse_he*=0.5
    #coarse_he*=0.5
    #coarse_he*=0.5
    coarse_area = coarse_he**2 / 2.0
    #fine_he = 0.5*coarse_he
    fine_he=coarse_he
    he = fine_he
    fine_area = fine_he**2 / 2.0
    boundaries=['left','right','bottom','top','front','back','obstacle']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    if structured:
        nnx=ceil(L[0]/he)+1
        nny=ceil(L[1]/he)+1
    else:
        vertices=[[maxX,inflowHeightMean-5.0*amplitude],#0
                  [maxX,inflowHeightMean+5.0*amplitude],#1
                  [maxX,minZ+Lz],#2
                  [transect[benchEnd][0],minZ+Lz],#3
                  [transect[benchStart][0],minZ+Lz],#4
                  [transect[0,0],minZ+Lz],#5
                  [transect[0,0],inflowHeightMean+5.0*amplitude],#6
                  [transect[0,0],inflowHeightMean-5.0*amplitude]]#7
        vertexFlags=[boundaryTags['right'],#0
                     boundaryTags['right'],#1
                     boundaryTags['top'],#2
                     boundaryTags['top'],#3
                     boundaryTags['top'],#4
                     boundaryTags['top'],#5
                     boundaryTags['right'],#6
                     boundaryTags['right']]#7
        for p in transect:
            vertices.append([p[0],p[1]])
            vertexFlags.append(boundaryTags['bottom'])
        nSegments = len(vertices)
        segments=[]
        segmentFlags=[]
        for i in range(nSegments):
            segments.append([i,(i+1) % nSegments])
            if i in [0,1,nSegments-1]:
                segmentFlags.append(boundaryTags['right'])
            elif i in [2,3,4]:
                segmentFlags.append(boundaryTags['top'])
            elif i in [5,6,7]:
                segmentFlags.append(boundaryTags['left'])
            else:
                segmentFlags.append(boundaryTags['bottom'])
        #mesh grading
        vertices.append([transect[benchStart][0],inflowHeightMean+5.0*amplitude])
        vertexFlags.append(0)
        ul=nSegments
        vertices.append([transect[benchStart][0],inflowHeightMean-5.0*amplitude])
        vertexFlags.append(0)
        ll=nSegments+1
        segments.append([ul,4])
        segmentFlags.append(0)
        segments.append([ul,6])
        segmentFlags.append(0)
        segments.append([ll,7])
        segmentFlags.append(0)
        segments.append([ll,benchStart+8])
        segmentFlags.append(0)
        eps=1.0e-5
        #regions=[[transect[0,0]+eps,transect[0,1]+eps],
        #         [transect[0,0]+eps,inflowHeightMean+5.0*amplitude+eps],
        #         [transect[0,0]+eps,inflowHeightMean-5.0*amplitude+eps]]
        #regionFlags=[1,2,3]
        #regionConstraints=[coarse_area,coarse_area,fine_area]
        domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                                      vertexFlags=vertexFlags,
                                                      segments=segments,
                                                      segmentFlags=segmentFlags)
        #regions=regions,
         #                                             regionFlags=regionFlags,
          #                                            regionConstraints=regionConstraints)
        #go ahead and add a boundary tags member 
        domain.boundaryTags = boundaryTags
        domain.writePoly("mesh")
        domain.writePLY("mesh")
        domain.writeAsymptote("mesh")
        triangleOptions="VApq30Dena%8.8f" % ((he**2)/2.0,)
        #triangleOptions="VApq30Dena"#%8.8f" % ((he**2)/2.0,)
# Time stepping
T=100*period
dt_fixed =period/10.0
dt_init = min(0.01*dt_fixed,0.001)
#runCFL=0.33
runCFL=0.033
nDTout = int(round(T/dt_fixed))

# Numerical parameters
ns_forceStrongDirichlet = False
if useMetrics:
    ns_shockCapturingFactor  = 0.9
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.9
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.5
    vof_shockCapturingFactor = 0.9
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.5
    rd_shockCapturingFactor  = 0.9
    rd_lag_shockCapturing = False
    epsFact_density    = 1.5
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 10.0
    redist_Newton = True
    kappa_shockCapturingFactor = 0.9
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

ns_nl_atol_res = max(1.0e-12,0.001*he**2)
vof_nl_atol_res = max(1.0e-12,0.001*he**2)
ls_nl_atol_res = max(1.0e-12,0.001*he**2)
rd_nl_atol_res = max(1.0e-12,0.1*he)
mcorr_nl_atol_res = max(1.0e-12,0.001*he**2)
kappa_nl_atol_res = max(1.0e-12,0.001*he**2)
dissipation_nl_atol_res = max(1.0e-12,0.001*he**2)

#turbulence
ns_closure=2 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
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

# Initial condition
waterLine_x = 215.15#2*L[0]
waterLine_z = inflowHeightMean#4m above bench
#waterLine_x = 0.5*L[0]
#waterLine_z = 0.9*L[1]

def signedDistance(x):
    phi_x = x[0]-waterLine_x
    phi_z = x[1]-waterLine_z 
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


def theta(x,t):
    return k*x[0] - omega*t

def z(x):
    return x[1] - inflowHeightMean

sigma = omega - k*inflowVelocityMean[0]
h = inflowHeightMean - transect[0][1]

def waveHeight(x,t):
    return inflowHeightMean + amplitude*cos(theta(x,t))

def waveVelocity_u(x,t):
    return sigma*amplitude*cosh(k*(z(x)+h))*cos(theta(x,t))/sinh(k*h)

def waveVelocity_v(x,t):
    return sigma*amplitude*sinh(k*(z(x)+h))*sin(theta(x,t))/sinh(k*h)

#solution variables

def wavePhi(x,t):
    return x[1] - waveHeight(x,t)

def waveVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t))

def twpflowVelocity_u(x,t):
    waterspeed = waveVelocity_u(x,t)
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
    u = H*windVelocity[0] + (1.0-H)*waterspeed
    return u

def twpflowVelocity_v(x,t):
    waterspeed = 0.0
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
    return H*windVelocity[1]+(1.0-H)*waterspeed

def twpflowFlux(x,t):
    return -twpflowVelocity_u(x,t)

def outflowVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,x[1] - outflowHeight)
