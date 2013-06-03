from math import *
import proteus.MeshTools
from proteus import Domain
from proteus import AnalyticalSolutions
from proteus.default_n import *   
   
#----------------------------------------------------
#  Discretization -- input options    
#----------------------------------------------------
Refinement=4
genMesh=True
spaceOrder=1
useRBLES   = 0.0
useMetrics = 1.0
use_petsc4py=False
useVF=1.0
#type of 2 equation turbulence model to use
# 0 None
#1 K-Epsilon
#2 Wilcox K-Omega, 1998
#3 Wilcox K-Omega, 1988
useRANS = 0 
            
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
    basis=C0_AffineLinearOnSimplexWithNodalBasis
    elementQuadrature = SimplexGaussQuadrature(nd,2)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,2) 	    
elif spaceOrder == 2:
    hFactor=0.5
    basis=C0_AffineQuadraticOnSimplexWithNodalBasis	
    elementQuadrature = SimplexGaussQuadrature(nd,4)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)


#----------------------------------------------------
# Domain and mesh
# 
#  
# ------------------------------------- L_z
# | ---> v^a_{in}                     | 
# |                                   |    p= p_{out}
# ************************************* -- z_g
# |          \theta_s                 | 
# |                                   |
# ------------------------------------- L_x
#
#----------------------------------------------------
#----------------------------------------------------
# Domain and mesh
# 
#               top                   L_z
# ------------------------------------- L_x
# | ---> v^a_{in}                     | 
# |                                   |    p= p_{out}
# ************************************* -- z_g
# l_duct|        porous          |r_duct
#       |                        |
#       ------------------------- L_px
#      0,0       bottom
#----------------------------------------------------

# Domain and mesh [m]
left_duct_length= 2.0
right_duct_length= 2.0
porous_length = 6.0   #sand box length
duct_length = left_duct_length+porous_length+right_duct_length  #free-flow domain length
z_g = 1.0 #height of tank
duct_height = 1.0
upstream_height = z_g+duct_height #total height

#bounding box
length = duct_length  
L = (length,
     upstream_height)


he = L[1]/float(4*Refinement-1)

#--- define Piecewise Linear Complex Domain --- 
boundaries=['left_fluid','right_fluid','bottom_left','bottom','bottom_right','top','left_porous','right_porous','interior_porous']
boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])

vertices=[[0.0,0.0], #0
          [porous_length,0.0],#1
          [porous_length,z_g],#2
          [porous_length+right_duct_length,z_g], #3
          [porous_length+right_duct_length,z_g+duct_height],#4
          [-left_duct_length,z_g+duct_height], #5
          [-left_duct_length,z_g],  #6
          [0.0,z_g]] #7

vertexFlags=[boundaryTags['bottom'],
             boundaryTags['bottom'],
             boundaryTags['right_porous'],
             boundaryTags['bottom'],
             boundaryTags['right_fluid'],
             boundaryTags['top'],
             boundaryTags['left_fluid'],
             boundaryTags['bottom']]

segments=[[0,1],
          [1,2],
          [2,3],
          [3,4],
          [4,5],
          [5,6],
          [6,7],
          [7,0],
          [2,7]]

segmentFlags=[boundaryTags['bottom'],
              boundaryTags['right_porous'],
              boundaryTags['bottom_right'],
              boundaryTags['right_fluid'],
              boundaryTags['top'],
              boundaryTags['left_fluid'],
              boundaryTags['bottom_left'],
              boundaryTags['left_porous'],
              boundaryTags['interior_porous']]

fluid_id = 0; porous_id = 1
regions=[[0.+0.5*porous_length,z_g+0.5*duct_height],[0.+0.5*porous_length,0.5*z_g]]
regionFlags=[fluid_id,porous_id]
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

#homogeneous subsurface
meanGrainSize = 0.001
porosity  = 0.3
porosityTypes = numpy.array([1.0,porosity])
#convert Ergun formula to alpha-beta coefficients
if porosity > 9.999e-1:
    Ftilde = 1.0
    Kinv   = 0.0
else:
    Ftilde = porosity*meanGrainSize*100.0/(1.0-porosity) #viscocities cancel out
    Kinv   = 180.*(1.0-porosity)*(1.0-porosity)/(meanGrainSize*meanGrainSize*porosity*porosity*porosity);
dragAlphaTypes = numpy.array([0.0,Kinv])
dragBetaTypes  = numpy.array([0.0,Kinv*Ftilde])

killNonlinearDragInSpongeLayer = False 

#--- Numerical parameters for mesh
nLevels = 1
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0

#----------------------------------------------------
# Physical coefficients
#----------------------------------------------------
# Air
rho_1 = 1.205
nu_1  = 1.500e-5 

#single phase
rho_0 = rho_1
nu_0  = rho_1

nu = nu_1
rho = rho_1

#Re = 10000.0#30250.0
#inflow = nu*Re/upstream_height #perhaps should be nu*Re/(upstream_height-z_g)

inflow = 1.0 
Re = inflow*upstream_height/nu

# Surface tension
sigma_01 = 0.0

# Gravity
g = [0.0,-9.8]

#turbulence
ns_closure=0 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
if useRANS == 1:
    ns_closure = 3
elif useRANS == 2:
    ns_closure == 4

#----------------------------------------------------
# Time stepping 
#----------------------------------------------------

residence_time = length/inflow
T=10.0*residence_time
tnList = [0.0,0.0001]
nDTout=10#100
tnList.extend([max(i*T/float(nDTout),0.1) for i in range(1,nDTout+1)])#[0.0,0.5*T,T]


#----------------------------------------------------
# Boundary conditions and other flags
#----------------------------------------------------
weak_bc_penalty_constant = 100.0


upstream_start_z = z_g
kInflow = 0.03*inflow*inflow#0.003*inflow*inflow
inflowVelocity = (inflow,0.0)

class u_flat:
    def __init__(self,val=inflow,ztop=upstream_height,zbot=upstream_start_z,delta_z=0.1):
        self.val= val; self.ztop = ztop; self.zbot=zbot
        self.delta_z = delta_z
    def uOfX(self,x):
        fact = 0.0#exp(-(x[1]-self.zbot)*(self.ztop-x[1])/self.delta_z)
        return self.val*(1.0-fact)

uProfile = u_flat(val=inflow)

def velRamp(t):
    if t < residence_time:
        return 1.0-exp(-25.0*t/residence_time)
    else:
        return 1.0

def twpflowPressure(x,flag):
    if flag in [boundaryTags['right_fluid'],boundaryTags['right_porous']]:#set pressure on outflow to hydrostatic
        return lambda x,t: -(L[1]-x[1])*rho*g[1]


bottom = [boundaryTags['bottom'],boundaryTags['bottom_left'],boundaryTags['bottom_right']]
porous_boundary = [boundaryTags['left_porous'],boundaryTags['right_porous']]
exterior_boundary = [boundaryTags['left_fluid'],boundaryTags['right_fluid'],boundaryTags['top'],boundaryTags['left_porous'],boundaryTags['right_porous']] + bottom
def twpflowVelocity_u(x,flag):
    if flag == boundaryTags['left_fluid']:
        return lambda x,t: inflow*velRamp(t)
    elif (flag == boundaryTags['top'] or
          flag in bottom or
          flag == boundaryTags['left_porous'] or 
          flag == boundaryTags['right_porous']):
        return lambda x,t: 0.0
def twpflowVelocity_v(x,flag):
    if flag in [boundaryTags['left_fluid'],boundaryTags['left_porous'],boundaryTags['right_porous']]:
        return lambda x,t: 0.0
    if (flag == boundaryTags['top'] or
        flag in bottom):
        return lambda x,t: 0.0


#----------------------------------------------------
# Initial condition
#----------------------------------------------------
#free surface
def signedDistance(x):
    phi_z = x[1]-L[1]*10.0 # (single phase )
    return phi_z

#----------------------------------------------------
# Numerical parameters
#----------------------------------------------------
# Numerical parameters
ns_forceStrongDirichlet = True
epsFact_solid = 0.5
if useMetrics:
    ns_shockCapturingFactor  = 0.1
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.1
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.0
    vof_shockCapturingFactor = 0.1
    vof_lag_shockCapturing = True
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
    kappa_lag_shockCapturing = True
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.1
    dissipation_lag_shockCapturing = True
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.0
else:
    ns_shockCapturingFactor  = 0.9
    ns_lag_shockCapturing = False
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.9
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.0
    rd_shockCapturingFactor  = 0.9
    rd_lag_shockCapturing = False
    epsFact_density    = 1.5
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 10.0
    redist_Newton = False
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True
    kappa_sc_uref  = 1.0
    kappa_sc_beta  = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True
    dissipation_sc_uref  = 1.0
    dissipation_sc_beta  = 1.0

ns_nl_atol_res = max(1.0e-8,0.1*he**2/2.0)
ls_nl_atol_res = max(1.0e-8,0.1*he**2/2.0)
rd_nl_atol_res = max(1.0e-8,0.1*he)
kappa_nl_atol_res = max(1.0e-8,0.01*he**2/2.0)
dissipation_nl_atol_res = max(1.0e-8,0.01*he**2/2.0)

#================================================
# Print run data
#================================================
logEvent("""
Reynolds number    = %16.21e
Windtunnel Domain  = [%16.21e,%16.21e]
Porous Domain      = [%16.21e,%16.21e]
Target element h_e = %16.21e
Porosity           = %16.21e
MeanGrainSize      = %16.21e
Inflow velocity    = %16.21e
T                  = %16.21e
nDTout             = %i
""" % (Re,           
       L[0],L[1],
       L[0],z_g,
       he,
       porosity,
       meanGrainSize,
       inflow,
       T,
       nDTout))
