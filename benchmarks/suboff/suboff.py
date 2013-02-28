from math import *
import proteus.MeshTools
from proteus import Domain
from proteus import AnalyticalSolutions
from proteus.default_n import *   
import suboff_domain
   
#----------------------------------------------------
#  Discretization -- input options    
#----------------------------------------------------
Refinement=1
genMesh=True
spaceOrder=1
useHex=False
useRBLES   = 0.0
useMetrics = 0.0
use_petsc4py=True
quasi2D = True
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


#----------------------------------------------------
# Domain and mesh
#----------------------------------------------------

if useHex: 
    assert False, "useHex not enabled yet for suboff problem"
    hex=True 

    comm=Comm.get()	
    if comm.isMaster():	
        size = numpy.array([[0.520,0.510   ,0.520],
	                    [0.330,0.335833,0.330],
			    [0.320,0.325   ,0.000]])/float(Refinement)
        numpy.savetxt('size.mesh', size)
        failed = os.system("../../scripts/marinHexMesh")      
     
    domain = Domain.MeshHexDomain("marinHex") 
else:
    #discretization size for theta in y,z plane
    ntheta = 8
    #max number of points in x, need to get rid of this
    max_nx = 300
    #how much to pad bounding box
    pad_x = 8.; pad_r_fact = 6.0
    #convert from feet to meters
    feet2meter = 0.3048
    #returns x,y,z points along hull, 
    #number of poins along y,z plane for each x
    #total number of points used in x
    #start of bounding box
    #size of bounding box
    x,y,z,ntheta_user,np,x_ll,L = suboff_domain.darpa2gen(max_nx,ntheta,
                                                          pad_x=pad_x,pad_r_fact=pad_r_fact,
                                                          length_conv=feet2meter)
    #just the points
    suboff_domain.write_csv_file(x[:np],y[:np],z[:np],'darpa2')
    #get the domain
    domain = suboff_domain.build_domain_from_axisymmetric_points(x[:np,:],y[:np,:],z[:np,:],x_ll,L,include_front_and_back=0,ntheta_user=ntheta_user[:np],name='darpa2')
    boundaryTags = domain.boundaryTags
    
    upstream_height=L[2]
    width = L[1]
    length = L[0]
    sub_length = L[0]-pad_x
    he = (sub_length)/float(6.5*Refinement)
    if quasi2D:
        width = he

    domain.writePLY('mesh_darpa2')
    domain.writePoly('mesh_darpa2')
    domain.writeAsymptote("mesh")
    triangleOptions="VApq1.25q12ena%e" % ((he**3)/6.0,)
#
nLevels = 1
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
nLayersOfOverlapForParallel = 0

#----------------------------------------------------
# Physical coefficients
#----------------------------------------------------
Re = 302.5#3025.0
inflow = 1.0
nu = inflow*sub_length/Re

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
g = [0.0,0.0,0.0]


#----------------------------------------------------
# Time stepping 
#----------------------------------------------------

residence_time = length/inflow
T=10.0#10.0*length/inflow
#tnList = [0.0,0.1*residence_time,T]
nDTout=10
tnList = [i*T/float(nDTout) for i in range(nDTout+1)]#[0.0,0.5*T,T]


#----------------------------------------------------
# Boundary conditions and other flags
#----------------------------------------------------
grad_p = -inflow/(upstream_height**2/(8.0*mu_0))
upstream_start_z = x_ll[2]
kInflow = 0.003*inflow*inflow

class u_flat:
    def __init__(self,val=inflow,ztop=upstream_height,zbot=0.0,delta_z=0.1):
        self.val= val; self.ztop = ztop; self.zbot=zbot
        self.delta_z = delta_z
    def uOfX(self,x):
        fact = exp(-(x[2]-self.zbot)*(self.ztop-x[2])/self.delta_z)
        return self.val*(1.0-fact)

uProfile = u_flat(val=inflow)

def velRamp(t):
    if t < residence_time:
        return 1.0-exp(-25.0*t/residence_time)
    else:
        return 1.0

def twpflowPressure(x,flag):
    if flag == boundaryTags['right']:#set pressure on outflow to hydrostatic
        return lambda x,t: -(downstream_height-x[2])*rho_0*g[2]


bottom = [boundaryTags['bottom']]

def twpflowVelocity_u(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: uProfile.uOfX(x-[0.0,0.0,upstream_start_z])*velRamp(t)
    elif (flag == boundaryTags['top'] or
          flag in bottom):
        return lambda x,t: 0.0

def twpflowVelocity_v(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: 0.0
    if (flag == boundaryTags['top'] or
        flag in bottom):
        return lambda x,t: 0.0

def twpflowVelocity_w(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: 0.0
    if (flag == boundaryTags['top'] or
        flag in bottom):
        return lambda x,t: 0.0

#----------------------------------------------------
# Initial condition
#----------------------------------------------------
waterLine_z = L[2]*10.0
def signedDistance(x):
    phi_z = x[2]-waterLine_z 
    return phi_z

#----------------------------------------------------
# Numerical parameters
#----------------------------------------------------

ns_shockCapturingFactor=0.3

ls_shockCapturingFactor=0.3
ls_sc_uref = 1.0
ls_sc_beta = 1.5

vof_shockCapturingFactor=0.3
vof_sc_uref = 1.0
vof_sc_beta = 1.5

rd_shockCapturingFactor=0.3

kappa_shockCapturingFactor=0.9
kappa_sc_uref = 1.0
kappa_sc_beta = 1.5

epsilon_shockCapturingFactor=0.9
epsilon_sc_uref = 1.0
epsilon_sc_beta = 1.5

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
epsFact_consrv_diffusion=10.0
epsFact_vof              = epsFact 


