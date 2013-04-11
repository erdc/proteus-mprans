from math import *
import proteus.MeshTools
from proteus import Domain
from proteus import AnalyticalSolutions
from proteus.default_n import *   
import suboff_domain
   
#----------------------------------------------------
#  Discretization -- input options    
#----------------------------------------------------
Refinement=4
genMesh=True
spaceOrder=1
useRBLES   = 0.0
useMetrics = 1.0
use_petsc4py=False
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
#----------------------------------------------------

#discretization size for theta in y,z plane
ntheta = 2
#max number of points in x, need to get rid of this
max_nx = 300
#how much to pad bounding box
pad_x = 10.; pad_r_fact = 4.0
#convert from feet to meters
feet2meter = 0.3048
#returns x,y,z points along hull, 
#number of poins along y,z plane for each x
#total number of points used in x
#start of bounding box
#size of bounding box
x,y,z,theta_offset,np,x_ll,L = suboff_domain.darpa2gen(max_nx,ntheta,
                                                       pad_x=pad_x,pad_r_fact=pad_r_fact,
                                                       length_conv=feet2meter)

L[2] = 1.0
#just the points
suboff_domain.write_csv_file(x[:np],y[:np],z[:np],'darpa2_2d')
#get the domain
domain = suboff_domain.build_2d_domain_from_axisymmetric_points(x[:np,:],y[:np,:],x_ll,L,name='darpa2_2d')
boundaryTags = domain.boundaryTags
    
upstream_height=L[1]
length = L[0]
sub_length = L[0]-pad_x
he = (sub_length)/float(6.5*Refinement)

domain.writePoly('mesh_darpa2_2d')
triangleOptions="VApq30Dena%8.8f" % ((he**2)/2.0,)
#
nLevels = 1
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
nLayersOfOverlapForParallel = 0

#----------------------------------------------------
# Physical coefficients
#----------------------------------------------------
Re = 302500.0
nu_h20 = 1.004e-6
nu = nu_h20
inflow = nu*Re/upstream_height


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
g = [0.0,-9.8]


#----------------------------------------------------
# Time stepping 
#----------------------------------------------------

residence_time = length/inflow
T=1.0*residence_time
tnList = [0.0,0.1]
nDTout=10
tnList.extend([max(i*T/float(nDTout),0.1) for i in range(1,nDTout+1)])#[0.0,0.5*T,T]


#----------------------------------------------------
# Boundary conditions and other flags
#----------------------------------------------------
grad_p = -inflow/(upstream_height**2/(8.0*mu_0))
upstream_start_z = x_ll[1]
kInflow = 0.003*inflow*inflow
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
    if flag == boundaryTags['right']:#set pressure on outflow to hydrostatic
        return lambda x,t: -(L[1]-x[1])*rho_0*g[1]


bottom = [boundaryTags['bottom']]

def twpflowVelocity_u(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: uProfile.uOfX(x-[0.0,upstream_start_z,0.0])*velRamp(t)
    elif (flag == boundaryTags['top'] or
          flag in bottom or
          flag == boundaryTags['obstacle']):
        return lambda x,t: 0.0
def twpflowVelocity_v(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: 0.0
    if (flag == boundaryTags['top'] or
        flag in bottom or
        flag == boundaryTags['obstacle']):
        return lambda x,t: 0.0


#----------------------------------------------------
# Initial condition
#----------------------------------------------------
waterLine_z = L[1]*10.0
def signedDistance(x):
    phi_z = x[1]-waterLine_z 
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

