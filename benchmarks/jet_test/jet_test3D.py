from math import *
import proteus.MeshTools
from proteus import Domain
from proteus import AnalyticalSolutions
from proteus.default_n import *   
   
#----------------------------------------------------
#  Discretization -- input options    
#----------------------------------------------------
Refinement=4
genMesh=False
spaceOrder=1
useRBLES   = 0.0
useMetrics = 1.0
use_petsc4py=False
#turbulence models
turbulenceClosureModel=1 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega

#type of 2 equation turbulence model to use
#1 K-Epsilon
#2 Wilcox K-Omega, 1998
#3 Wilcox K-Omega, 1988
dissipation_model_flag = 0
if turbulenceClosureModel > 2:
    useRANS=1
else:
    useRANS=0


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
#----------------------------------------------------
he = 0.025 #guess level 0
domain = Domain.MeshTetgenDomain(fileprefix="jet_domain3d")
boundaryTags = {'inflow':5,'top_walls':6,'outlet':7,'outer_walls':8,'bottom_walls':9}
wall_boundaries= [boundaryTags['top_walls'],boundaryTags['outer_walls'],boundaryTags['bottom_walls']]
domain.boundaryTags = boundaryTags

L = (14.1,30.,30)
x_ll = (0.,-15.,-15.)
length = L[0]
inflow_length = 4.0*2. #diameter 
inflow_start_z= 0.0


#
nLevels = 1
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
nLayersOfOverlapForParallel = 0

#----------------------------------------------------
# Physical coefficients
#----------------------------------------------------
inflow = -4.0e-2
nu_h20 = 1.004e-6
nu = nu_h20
Re = length*abs(inflow)/nu


# Water
rho_0 = 998.2
nu_0  = nu
mu_0 = rho_0*nu_0

# Air
rho_1 = 1.205
nu_1  = 1.500e-5 

# Surface tension
sigma_01 = 0.0

# Gravity
g = [0.0,0.0,0.0]#[-9.8,0.0,0.0]


#----------------------------------------------------
# Time stepping 
#----------------------------------------------------

T=1.0
residence_time=1.
tnList = [0.0,0.0001]
nDTout=100
tnList.extend([max(i*T/float(nDTout),0.1) for i in range(1,nDTout+1)])#[0.0,0.5*T,T]


#----------------------------------------------------
# Boundary conditions and other flags
#----------------------------------------------------
kInflow = 0.03*inflow*inflow#0.003*inflow*inflow
inflowVelocity = (inflow,0.0,0.)

class u_flat:
    def __init__(self,val=inflow,ztop=inflow_length,zbot=inflow_start_z,delta_z=0.1):
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
    if flag == boundaryTags['outlet']:#set pressure on outflow to atmospheric
        return lambda x,t: 0.0



def twpflowVelocity_u(x,flag):
    if flag == domain.boundaryTags['inflow']:
        return lambda x,t: uProfile.uOfX(x-[0.0,inflow_start_z,0.0])*velRamp(t)
    if flag in wall_boundaries:
        return lambda x,t: 0.0
def twpflowVelocity_v(x,flag):
    if flag == boundaryTags['inflow']:
        return lambda x,t: 0.0
    if flag in wall_boundaries:
        return lambda x,t: 0.0
def twpflowVelocity_w(x,flag):
    if flag == boundaryTags['inflow']:
        return lambda x,t: 0.0
    if flag in wall_boundaries:
        return lambda x,t: 0.0


#----------------------------------------------------
# Initial condition
#----------------------------------------------------
waterLine_z = L[0]*10.0
def signedDistance(x):
    phi_z = x[0]-waterLine_z 
    return phi_z

#----------------------------------------------------
# Numerical parameters
#----------------------------------------------------
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
    ns_lag_shockCapturing = False
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.9
    ls_lag_shockCapturing = True#False
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
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref  = 1.0
    kappa_sc_beta  = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref  = 1.0
    dissipation_sc_beta  = 1.0

ns_nl_atol_res = max(1.0e-8,0.1*he**2/2.0)
ls_nl_atol_res = max(1.0e-8,0.1*he**2/2.0)
rd_nl_atol_res = max(1.0e-8,0.1*he)
kappa_nl_atol_res = max(1.0e-8,0.01*he**2/2.0)
dissipation_nl_atol_res = max(1.0e-8,0.01*he**2/2.0)

