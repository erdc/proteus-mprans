from math import *
import proteus.MeshTools

checkMass=False
applyCorrection=True
applyRedistancing=True
rdtimeIntegration='osher'
freezeLevelSet=True#False

spaceOrder=1
useHex=False
if spaceOrder == 1:
    sloshbox_quad_order = 2
elif spaceOrder == 2:
    sloshbox_quad_order = 4

nodalPartitioning=True
nd = 3
useBackwardEuler=True
useBackwardEuler_ls=True
useFixedStep=True
runCFL = 0.33
dt_fixed = 0.01
timeOrder = 1
dt_init=0.01
T=0.40
nDTout = int(round(T/dt_fixed))

lag_ns_subgridError=False
lag_ns_shockCapturing=False
lag_ls_shockCapturing=False

ns_shockCapturingFactor=0.2
ls_shockCapturingFactor=0.2
vof_shockCapturingFactor=0.2
rd_shockCapturingFactor=0.2

#epsilons for Heaviside/Dirac/etc smoothing
epsFact_density = 1.5
epsFact_viscosity = 1.5
epsFact_redistance = 0.33
epsFact_curvature=1.5
epsFact_consrv_heaviside=1.5
epsFact_consrv_dirac=1.5
epsFact_consrv_diffusion=10.0
epsFact_vof=1.5

usePETSc=True

L = (0.584,0.146,0.350)
useStokes=False

nLevels = 1#
useShock=True

lRefinement=2
nnx=(2**lRefinement)*10+1
nnx=40
nny=10
nnz=20
he = L[0]/float(nnx-1)

if useHex:
    hex=True
    soname="sloshbox_c0q"+`spaceOrder`+"_bdf_"+`timeOrder`+"_level_"+`lRefinement`
else:
    soname="sloshbox_c0p"+`spaceOrder`+"_bdf_"+`timeOrder`+"_level_"+`lRefinement`

pseudo2D=False
if pseudo2D:
    nny=2
    L=(L[0],he,L[2])
    shock_x = 0.5*L[0]
    shock_y = -L[0]#hack to get signed distance right in 2D
    shock_z = 0.5*L[2]
else:
    #nny = int((L[1]/he)+0.5) + 1
    #nnz = int((L[2]/he)+0.5) + 1
    shock_x = L[0]-0.146
    shock_y = -20.0   #L[1]
    shock_z = 0.292

print he
print nnx,nny,nnz
print  L[0]/float(nnx-1), L[1]/float(nny-1), L[2]/float(nnz-1)

regularGrid=True
if regularGrid:
    from proteus.Domain import RectangularDomain
    domain = RectangularDomain(L)
else:
    from tank3dDomain import *
    domain = tank3d(L=L)
    domain.writePoly("tank3d")
    domain.writePLY("tank3d")
    domain.writeAsymptote("tank3d")
    triangleOptions="VApq2q10ena%21.16e" % ((he**3)/6.0,)
    bt = domain.boundaryTags
    
waterLevel = 0.292
slopeAngle = 0.5*(pi/2.0)

#water
rho_0=998.2
nu_0=1.004e-6
#air
rho_1=1.205
nu_1= 1.500e-5 #* 1000.0
#rho_1 = rho_0
#nu_1 = nu_0
sigma_01=0.0#72.8e-3
#gravity
g=[0.0,0.0,-9.8]
closedTop=False
sidesNoSlip=False
bcsTimeDependent = False
def getDBC_p_sloshbox(x,flag):
    if closedTop:        
        if (x[0] < 1.0e-8 and x[2] > L[2] - 1.0e-8):#top left corner
            return lambda x,t: 0.0
    else:
        if x[2] > L[2] - 1.0e-8:
                return lambda x,t: 0.0


def getDBC_u_sloshbox(x,flag):
    return None

def getDBC_v_sloshbox(x,flag):
    return None

def getDBC_w_sloshbox(x,flag):
    return None

dirichletConditions = {0:getDBC_p_sloshbox,
                       1:getDBC_u_sloshbox,
                       2:getDBC_v_sloshbox,
                       3:getDBC_w_sloshbox}

def getAFBC_p_sloshbox(x,flag):
    if closedTop:
            if(x[0] == 0.0 or
               x[0] == L[0] or
               x[1] == 0.0 or
               x[1] == L[1] or
               x[2] == 0.0 or
               x[2] == L[2]):
                return lambda x,t: 0.0
    else:
            if x[2] > L[2] - 1.0e-8:
                return None
            else:
                return lambda x,t: 0.0

def getAFBC_u_sloshbox(x,flag):
    if (x[0] == 0.0 or
        x[0] == L[0]):
       return lambda x,t: 0.0

def getAFBC_v_sloshbox(x,flag):
    if (x[1] == 0.0 or
        x[1] == L[1]):
       return lambda x,t: 0.0

def getAFBC_w_sloshbox(x,flag):
    if  (x[2] == 0.0 or
        (x[2] == L[2] and closedTop)  ):
       return lambda x,t: 0.0

def getDFBC_u_sloshbox(x,flag):
    return lambda x,t: 0.0
    
def getDFBC_v_sloshbox(x,flag):
    return lambda x,t: 0.0

def getDFBC_w_sloshbox(x,flag):
    return lambda x,t: 0.0

fluxBoundaryConditions = {0:'mixedFlow',
                          1:'mixedFlow',
                          2:'mixedFlow',
                          3:'mixedFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_sloshbox,
                                    1:getAFBC_u_sloshbox,
                                    2:getAFBC_v_sloshbox,
                                    3:getAFBC_w_sloshbox}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u_sloshbox},
                                   2:{2:getDFBC_v_sloshbox},
                                   3:{3:getDFBC_w_sloshbox}}

def shockSignedDistance(x):
    phi_x = shock_x - x[0]
    phi_y = shock_y - x[1]
    phi_z = x[2] - shock_z
    phi_xy = sqrt((x[0] - shock_x)**2 + (x[1] - shock_y)**2)
    phi_xz = sqrt((x[0] - shock_x)**2 + (x[2] - shock_z)**2)
    phi_yz = sqrt((x[1] - shock_y)**2 + (x[2] - shock_z)**2)
    phi_xyz = sqrt((x[0] - shock_x)**2 + (x[1] - shock_y)**2 + (x[2] - shock_z)**2)
    if x[0] < shock_x:
        if x[2] < shock_z:
            if x[1] < shock_y:
                return phi_xy
            else:
                return phi_x
        else:
            if x[1] < shock_y:
                return phi_xyz
            else:
                return phi_xz
    elif x[1] < shock_y:
        if x[2] < shock_z:
            return phi_y
        else:
            return phi_yz
    else:
        if x[2] > shock_z:
            return phi_z
        else:
            return -min(fabs(phi_x),min(fabs(phi_y),fabs(phi_z)))

class PerturbedSurface_p:
    def __init__(self,waterLevel,slopeAngle):
        self.waterLevel=waterLevel
        self.slopeAngle=slopeAngle
    def uOfXT(self,x,t):
        if useShock:
            if shockSignedDistance(x) < 0:
                return -(L[2] - shock_z)*rho_1*g[2] - (shock_z - x[2])*rho_0*g[2]
            else:
                return -(L[2] - shock_z)*rho_1*g[2]
        else:
            surfaceNormal = [-sin(self.slopeAngle),cos(self.slopeAngle)]
            topRight_z= min(L[2],-(L[0]-0.5)*surfaceNormal[0] + self.waterLevel*surfaceNormal[1])
            signedDistance = (x[0] - 0.5)*surfaceNormal[0]+(x[2] - self.waterLevel)*surfaceNormal[1]
            signedDistanceCorner = (0.0 - 0.5)*surfaceNormal[0]+(L[2] - self.waterLevel)*surfaceNormal[1]
            if signedDistance > 0.0:
                p = -(signedDistanceCorner-signedDistance)*rho_1*g[2]
            else:
                p = -signedDistanceCorner*rho_1*g[2] + signedDistance*rho_0*g[2]
            return p

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:PerturbedSurface_p(waterLevel,slopeAngle),
                     1:AtRest(),
                     2:AtRest(),
                     3:AtRest()}

restrictFineSolutionToAllMeshes=False
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 1
#force true residual test
#linearSolverConvergenceTest = 'r-true' #r,its,r-true for true residual
