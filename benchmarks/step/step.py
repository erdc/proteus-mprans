from math import *
import proteus.MeshTools
from proteus import Domain
from proteus import AnalyticalSolutions
from proteus.default_n import *   
import step3d
   
#----------------------------------------------------
#  Discretization -- input options    
#----------------------------------------------------
Refinement=16
genMesh=True
spaceOrder=1
useHex=False
useRBLES   = 0.0
useMetrics = 1.0
use_petsc4py=True
quasi2D = True
use_PlanePoiseuille = False
use_KOmega = True
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
    assert False, "useHex not enabled yet for step problem"
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
    upstream_height=0.5
    width = 0.3*upstream_height
    downstream_height=1.0
    upstream_length = 1.0
    downstream_length = 5.0
    length = upstream_length+downstream_length
    polyfile = "step3d"
    he = length/float(6.5*Refinement)
    if quasi2D:
        width = he
    L      = [length,width,downstream_height]
    
    polyfile = "step3d"
    boundaryTags = step3d.genPoly(fileprefix=polyfile,
                                  width=width,
                                  upstream_height=upstream_height,
                                  downstream_height=downstream_height,
                                  upstream_length=upstream_length,
                                  downstream_length=downstream_length,
                                  step_fun=step3d.linear_profile,
                                  n_points_on_step=2,
                                  step_length=0.0)

    domain = Domain.PiecewiseLinearComplexDomain(fileprefix=polyfile)
    #go ahead and add a boundary tags member 
    domain.boundaryTags = boundaryTags
    domain.writePoly("mesh")
    domain.writePLY("mesh")
    domain.writeAsymptote("mesh")
    triangleOptions="VApq1.25q12ena%e" % ((he**3)/6.0,)
#
nLevels = 1
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
nLayersOfOverlapForParallel = 0

#----------------------------------------------------
# Physical coefficients
#----------------------------------------------------
Re = 20000.0#3025.0
nu_h20 = 1.004e-6
#inflow = 1.0
#nu = inflow*(downstream_height-upstream_height)/Re
nu = nu_h20
inflow = nu*Re/(downstream_height-upstream_height)

comm=Comm.get()
if comm.isMaster():
    print "step problem use_KOmega= %s Re= %12.5e nu =%12.5e " % (use_KOmega,Re,nu)
#
# Water
rho_0 = 998.2
nu_0  = 1.004e-6
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
T=60.0#10.0*length/inflow
#tnList = [0.0,0.1*residence_time,T]
nDTout=60
tnList = [0.0,0.00001]
#tnList = [i*T/float(nDTout) for i in range(nDTout+1)]#[0.0,0.5*T,T]
tnList.extend([i*T/float(nDTout) for i in range(1,nDTout+1)])#[0.0,0.5*T,T]
runCFL=0.9

#----------------------------------------------------
# Boundary conditions and other flags
#----------------------------------------------------
grad_p = -inflow/(upstream_height**2/(8.0*mu_0))
upstream_start_z = downstream_height - upstream_height
kInflow = 0.003*inflow*inflow

class u_flat:
    def __init__(self,val=inflow,ztop=upstream_height,zbot=0.0,delta_z=0.1):
        self.val= val; self.ztop = ztop; self.zbot=zbot
        self.delta_z = delta_z
    def uOfX(self,x):
        fact = exp(-(x[2]-self.zbot)*(self.ztop-x[2])/self.delta_z)
        return self.val*(1.0-fact)

if use_PlanePoiseuille:
    uProfile = AnalyticalSolutions.PlanePoiseuilleFlow_u2(plane_theta=0.0,
                                                          plane_phi=0.0,
                                                          v_theta=0.0,
                                                          v_phi=math.pi/2.0,
                                                          v_norm=0.0,
                                                          mu=mu_0,
                                                          grad_p=grad_p,
                                                          L=[1.0,1.0,upstream_height])

else:
    uProfile = u_flat(val=inflow)

def velRamp(t):
    if t < residence_time:
        return 1.0-exp(-25.0*t/residence_time)
    else:
        return 1.0

def twpflowPressure(x,flag):
    if flag == boundaryTags['downstream']:#set pressure on outflow to hydrostatic
        return lambda x,t: -(downstream_height-x[2])*rho_0*g[2]


bottom = [boundaryTags['upstream_bottom'],boundaryTags['step_bottom'],boundaryTags['downstream_bottom']]

def twpflowVelocity_u(x,flag):
    if flag == boundaryTags['upstream']:
        return lambda x,t: uProfile.uOfX(x-[0.0,0.0,upstream_start_z])*velRamp(t)
    elif (flag == boundaryTags['top'] or
          flag in bottom):
        return lambda x,t: 0.0

def twpflowVelocity_v(x,flag):
    if flag == boundaryTags['upstream']:
        return lambda x,t: 0.0
    if (flag == boundaryTags['top'] or
        flag in bottom):
        return lambda x,t: 0.0

def twpflowVelocity_w(x,flag):
    if flag == boundaryTags['upstream']:
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
ns_lag_shockCapturing = True
ns_lag_subgridError = True
if useMetrics:
    ns_shockCapturingFactor = 0.1
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

dissipation_shockCapturingFactor=0.9
dissipation_sc_uref = 1.0
dissipation_sc_beta = 1.5

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


