from math import *
import proteus.MeshTools
import numpy as np
#import waveModules_Matt as wm
from proteus import Domain
from proteus.default_n import *   
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
   
#  Discretization -- input options  
Refinement = 4#4#15
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
L = (5.0,
     10.0,
     1.0)
spongeLayer = True
xSponge = 8.0#0.8*L[0]
levee=True; spongeLayer=False;
leveeStart = 1.9
leveeBottomWidth = 3.0
leveeHeight =L[2]*3.0/5.0
leveeHeightDownstream = 0.0#0.25 
leveeSlope = 1.0/2.0
bedHeight = 0.2*L[2]
leveeHeightDownstream=bedHeight
quasi2D = True
veg=True; levee=False; spongeLayer=False
nLevels = 1
#parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0
rightEndClosed = False
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
    if spongeLayer:
        vertices=[[0.0,0.0,0.0],#0
                  [xSponge,0.0,0.0],#1
                  [xSponge,L[1],0.0],#2
                  [0.0,L[1],0.0],#3
                  [0.0,0.0,L[2]],#4
                  [xSponge,0.0,L[2]],#5
                  [xSponge,L[1],L[2]],#6
                  [0.0,L[1],L[2]],#7
                  [L[0],0.0,0.0],#8
                  [L[0],L[1],0.0],#9
                  [L[0],0.0,L[2]],#10
                  [L[0],L[1],L[2]]#11
                  ]
        vertexFlags=[boundaryTags['left'],
                     boundaryTags['front'],
                     boundaryTags['back'],
                     boundaryTags['left'],
                     boundaryTags['left'],
                     boundaryTags['front'],
                     boundaryTags['back'],
                     boundaryTags['left'],
                     boundaryTags['right'],
                     boundaryTags['right'],
                     boundaryTags['right'],
                     boundaryTags['right']]
        facets=[[[0,1,2,3]],#bottom
                [[0,1,5,4]],#front
                [[1,2,6,5]],#internal
                [[2,3,7,6]],#back
                [[3,0,4,7]],#left
                [[4,5,6,7]],#top
                [[1,8,9,2]],#bottom #start sponge
                [[1,8,10,5]],#front
                [[8,9,11,10]],#right
                [[2,6,11,9]],#back
                [[5,6,11,10]],#top
                ]
        facetFlags=[boundaryTags['bottom'],
                    boundaryTags['front'],
                    0,#boundaryTags['right'],
                    boundaryTags['back'],
                    boundaryTags['left'],
                    boundaryTags['top'],
                    boundaryTags['bottom'],#sponge
                    boundaryTags['front'],
                    boundaryTags['right'],
                    boundaryTags['back'],
                    boundaryTags['top']]
        regions=[[0.5*xSponge,0.5*L[1],0.5*L[2]],[0.5*(xSponge+L[0]),0.5*L[1],0.5*L[2]]]
        regionFlags=[0,1]
        spongeGrainSize= 0.05
        spongePorosity = 0.9
        killNonlinearDragInSpongeLayer = True#True
        porosityTypes      = numpy.array([1.0,spongePorosity])
        meanGrainSizeTypes = numpy.array([1.0,spongeGrainSize])
    elif levee:
        vertices=[[0.0,0.0,bedHeight],#0
                  [L[0],0.0,leveeHeightDownstream],#1
                  [L[0],L[1],leveeHeightDownstream],#2
                  [0.0,L[1],bedHeight],#3
                  [0.0,0.0,L[2]],#4
                  [L[0],0.0,L[2]],#5
                  [L[0],L[1],L[2]],#6
                  [0.0,L[1],L[2]],#7
                  [leveeStart,0.0,bedHeight],#8
                  [leveeStart+leveeBottomWidth,0.0,leveeHeightDownstream],#9
                  [leveeStart,L[1],bedHeight],#10
                  [leveeStart+leveeBottomWidth,L[1],leveeHeightDownstream],#11
                  [leveeStart+leveeHeight/leveeSlope,0.0,leveeHeight],#12
                  [leveeStart+leveeBottomWidth-leveeHeight/leveeSlope,0.0,leveeHeight],#13
                  [leveeStart+leveeHeight/leveeSlope,L[1],leveeHeight],#14
                  [leveeStart+leveeBottomWidth-leveeHeight/leveeSlope,L[1],leveeHeight],#15
                  [0.0,0.0,0.0],#16
                  [L[0],0.0,0.0],#17
                  [L[0],L[1],0.0],#18
                  [0.0,L[1],0.0]#19
                  ]
        vertexFlags=[boundaryTags['left'],
                     boundaryTags['right'],
                     boundaryTags['right'],
                     boundaryTags['left'],
                     boundaryTags['left'],
                     boundaryTags['right'],
                     boundaryTags['right'],
                     boundaryTags['left'],
                     0,#boundaryTags['bottom'],
                     0,#boundaryTags['bottom'],
                     0,#boundaryTags['bottom'],
                     0,#boundaryTags['bottom'],
                     boundaryTags['left'],
                     boundaryTags['left'],
                     boundaryTags['right'],
                     boundaryTags['right'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom']
                     ]
        facets=[[[0,8,10,3]],#bed
                [[8,9,11,10]],#bed
                [[9,1,2,11]],#bed
                [[8,9,13,12]],#front levee face
                [[10,11,15,14]],#back levee face
                [[8,12,14,10]],#left levee face
                [[12,13,15,14]],#top levee face
                [[13,9,11,15]],#right levee face
                [[0, 8,12,13, 9,1,5,4]], #front facet
                [[3,10,14,15,11,2,6,7]], #back facet
                [[1,2,6,5]],#right
                [[3,0,4,7]],#left
                [[4,5,6,7]],#top
                [[16,17,18,19]],#bottom
                [[16,17,1,9,8,0]],#front 
                [[19,18,2,11,10,3]],#back
                [[16,19,3,0]],#left
                [[17,18,2,1]]#right
                ]
        facetFlags=[0,#boundaryTags['bottom'],
                    0,#boundaryTags['bottom'],
                    0,#boundaryTags['bottom'],
                    boundaryTags['front'],
                    boundaryTags['back'],
                    0,
                    0,
                    0,
                    boundaryTags['front'],
                    boundaryTags['back'],
                    boundaryTags['right'],
                    boundaryTags['left'],
                    boundaryTags['top'],
                    boundaryTags['bottom'],
                    boundaryTags['front'],
                    boundaryTags['back'],
                    boundaryTags['left'],
                    boundaryTags['right']
                    ]
        regions=[[0.001,0.001,bedHeight+0.001],[leveeStart+0.5*leveeBottomWidth,0.5*L[1],0.5*(bedHeight+leveeHeight)],[0.5*L[0],0.5*L[1],0.5*bedHeight]]
        regionFlags=[0,1,2]
        spongeGrainSize= 0.003
        spongePorosity = 0.2
        killNonlinearDragInSpongeLayer = True
        porosityTypes      = numpy.array([1.0,spongePorosity,spongePorosity])
        meanGrainSizeTypes = numpy.array([1.0,spongeGrainSize,spongeGrainSize])
    elif veg:
        dh = 0.05 #m
        he = 0.5*dh
        L = (60*dh,3*dh,L[2])
        random_height=False
        if random_height:
            rg = random.Random()
        iSkip=30
        iStop=45
        nBoxes_xy=[60,3]
        boxScale=[0.5,0.5,0.5]
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
#               [box_xy[0],box_xy[1],0.0],#8
#               [box_xy[0]+box_L[0],box_xy[1],0.0],#9
#               [box_xy[0]+box_L[0],box_xy[1]+box_L[1],0.0],#10
#               [box_xy[0],box_xy[1]+box_L[1],0.0],#11
#               [box_xy[0],box_xy[1],box_L[2]],#12
#               [box_xy[0]+box_L[0],box_xy[1],box_L[2]],#13
#               [box_xy[0]+box_L[0],box_xy[1]+box_L[1],box_L[2]],#14
#               [box_xy[0],box_xy[1]+box_L[1],box_L[2]]]#15
        vertexFlags=[boundaryTags['left'],
                     boundaryTags['right'],
                     boundaryTags['right'],
                     boundaryTags['left'],
                     boundaryTags['left'],
                     boundaryTags['right'],
                     boundaryTags['right'],
                     boundaryTags['left']]
#                  boundaryTags['obstacle'],
#                  boundaryTags['obstacle'],
#                  boundaryTags['obstacle'],
#                  boundaryTags['obstacle'],
#                  boundaryTags['obstacle'],
#                  boundaryTags['obstacle'],
#                  boundaryTags['obstacle'],
#                  boundaryTags['obstacle']]
        facets=[[[0,1,2,3]],
                [[0,1,5,4]],
                [[1,2,6,5]],
                [[2,3,7,6]],
                [[3,0,4,7]],
                [[4,5,6,7]]]
        facetHoles=[[],
                    [],
                    [],
                    [],
                    [],
                    []]
#             [[8,9,13,12]],
#             [[9,10,14,13]],
#             [[10,11,15,14]],
#             [[11,8,12,15]],
#             [[12,13,14,15]]]
        facetFlags=[boundaryTags['bottom'],
                    boundaryTags['front'],
                    boundaryTags['right'],
                    boundaryTags['back'],
                    boundaryTags['left'],
                    boundaryTags['top']]
#                 boundaryTags['obstacle'],
#                 boundaryTags['obstacle'],
#                 boundaryTags['obstacle'],
#                 boundaryTags['obstacle'],
#                 boundaryTags['obstacle']]
        holes=[]
        dx = L[0]/float(nBoxes_xy[0])
        dy = L[1]/float(nBoxes_xy[1])
        hbdx = 0.5*boxScale[0]*dx
        hbdy = 0.5*boxScale[1]*dy
        for i in range(nBoxes_xy[0]):
            for j in range(nBoxes_xy[1]):
                if (i == nBoxes_xy[0] - 1 and j%2 == 1) or (i < iSkip) or (i >= iStop):
                    continue
                center = [i*dx + 0.5*dx+(j%2)*0.5*dx,j*dy+0.5*dy]
                hz = boxScale[2]*L[2]
                if random_height:
                    hz *= rg.gauss(1.0,0.1)
                #holes.append([center[0],center[1],0.5*hz])
                n=len(vertices)
                vertices.append([center[0]-hbdx,center[1]-hbdy,0.0])
                vertices.append([center[0]-hbdx,center[1]+hbdy,0.0])
                vertices.append([center[0]+hbdx,center[1]+hbdy,0.0])
                vertices.append([center[0]+hbdx,center[1]-hbdy,0.0])
                vertices.append([center[0]-hbdx,center[1]-hbdy,hz])
                vertices.append([center[0]-hbdx,center[1]+hbdy,hz])
                vertices.append([center[0]+hbdx,center[1]+hbdy,hz])
                vertices.append([center[0]+hbdx,center[1]-hbdy,hz])
                vertexFlags.append(boundaryTags['obstacle'])
                vertexFlags.append(boundaryTags['obstacle'])
                vertexFlags.append(boundaryTags['obstacle'])
                vertexFlags.append(boundaryTags['obstacle'])
                vertexFlags.append(boundaryTags['obstacle'])
                vertexFlags.append(boundaryTags['obstacle'])
                vertexFlags.append(boundaryTags['obstacle'])
                vertexFlags.append(boundaryTags['obstacle'])
                facets[0].append([n+0,n+1,n+2,n+3])#bottom
                fh = (center[0],center[1],0.0)
                facetHoles[0].append(fh)
                facets.append([[n+0,n+1,n+5,n+4]])#left
                facetHoles.append([])
                facets.append([[n+1,n+2,n+6,n+5]])#back
                facetHoles.append([])
                facets.append([[n+2,n+3,n+7,n+6]])#right
                facetHoles.append([])
                facets.append([[n+3,n+0,n+4,n+7]])#front
                facetHoles.append([])
                facets.append([[n+4,n+5,n+6,n+7]])#top
                facetHoles.append([])
                facetFlags.append(boundaryTags['obstacle'])
                facetFlags.append(boundaryTags['obstacle'])
                facetFlags.append(boundaryTags['obstacle'])
                facetFlags.append(boundaryTags['obstacle'])
                facetFlags.append(boundaryTags['obstacle'])        
        regions=[[0.001,0.001,0.001]]
        regionFlags=[1]
    else:
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
                                                 facetHoles=facetHoles,
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
ns_shockCapturingFactor  = 0.3
ls_shockCapturingFactor  = 0.3
ls_sc_uref  = 1.0
ls_sc_beta  = 1.0
vof_shockCapturingFactor = 0.3
vof_sc_uref = 1.0
vof_sc_beta = 1.0
rd_shockCapturingFactor  = 0.3

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
inflowVelocityMean = (0.0,0.0,0.0)

waveLength = inflowHeightMean #
period = waveLength/sqrt((-g[2])*inflowHeightMean) #meters
omega = 2.0*pi/period
k=(2.0*pi/waveLength,0.0,0.0)
amplitude = 0.1*inflowHeightMean

# Wave Field Object
#waveField = wm.Linear2D(amplitude,omega,k,L[2],rho_0,rho_1)
#waveField = wm.Solitary(amplitude,omega,k,L[2],rho_0,rho_1)

def waveHeight(x,t):
    return inflowHeightMean + amplitude*sin(omega*t-k[0]*x[0])
#    return inflowHeightMean + waveField.height(x,t)

def waveVelocity_u(x,t):
    return inflowVelocityMean[0] + omega*amplitude*sin(omega*t - k[0]*x[0])/(k[0]*inflowHeightMean)
#    z = x[2] - inflowHeightMean
#    return inflowVelocityMean[0] + waveField.velocity_u(x,t)

def waveVelocity_w(x,t):
    z = x[2] - inflowHeightMean
    return inflowVelocityMean[2] + (z + inflowHeightMean)*omega*amplitude*cos(omega*t-k[0]*x[0])/inflowHeightMean
#    return inflowVelocityMean[2] + waveField.velocity_w(x,t)
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
#    z = x[2] - inflowHeightMean
#    return waveField.pressure(x,t,z)

def outflowPressure(x,t):
    p_L = L[2]*rho_1*g[2]
    phi_L = L[2] - outflowHeight
    phi = x[2] - outflowHeight
    return p_L -g[2]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi_L)
                                                          -smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi)))

# Time 
T=20.0#period*100
runCFL = 0.33
print "T",T
dt_fixed = period/10.0 
#dt_fixed = period/100.0
#dt_fixed = 6.0/1000.0
dt_init = 1.0e-3
nDTout = int(T/dt_fixed)
tnList = [i*dt_fixed for i in range(0,nDTout+1)] 
print tnList
