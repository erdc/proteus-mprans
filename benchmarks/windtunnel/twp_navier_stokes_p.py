from proteus import *
from proteus.default_p import *
from windtunnel import *
from proteus.mprans import RANS2P

LevelModelType = RANS2P.LevelModel

if useRANS >= 1:
    Closure_0_model = 5; Closure_1_model=6
    if useOnlyVF:
        Closure_0_model=2; Closure_1_model=3
else:
    Closure_0_model = None
    Closure_1_model = None

coefficients = RANS2P.Coefficients(epsFact=epsFact_viscosity,
                                   sigma=0.0,
                                   rho_0 = rho,
                                   nu_0 = nu,
                                   rho_1 = rho,
                                   nu_1 = nu,
                                   g=g,
                                   nd=nd,
                                   LS_model=1,
                                   Closure_0_model=Closure_0_model,
                                   Closure_1_model=Closure_1_model,
                                   epsFact_density=epsFact_density,
                                   stokes=False,
                                   useVF=useVF,
				   useRBLES=useRBLES,
				   useMetrics=useMetrics,
                                   eb_adjoint_sigma=1.0,
                                   eb_penalty_constant=weak_bc_penalty_constant,
                                   forceStrongDirichlet=ns_forceStrongDirichlet,
                                   turbulenceClosureModel=ns_closure,
                                   porosityTypes=porosityTypes,
                                   dragAlphaTypes=dragAlphaTypes,
                                   dragBetaTypes=dragBetaTypes,
                                   epsFact_solid = epsFact_solid)


bcsTimeDependent = True

periodic = False

getDBC_p = twpflowPressure
getDBC_u = twpflowVelocity_u
getDBC_v = twpflowVelocity_v

dirichletConditions = {0:getDBC_p,
                       1:getDBC_u,
                       2:getDBC_v}


fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}


def getAFBC_p(x,flag):
    if flag == boundaryTags['left_fluid']: #inflow
        return lambda x,t: -inflow*velRamp(t)
    #if flag in [boundaryTags['right_fluid']]:#outflow
    #    return None
    if flag in [boundaryTags['top'],boundaryTags['bottom'],boundaryTags['left_porous'],
                boundaryTags['right_porous']]:
        return lambda x,t: 0.0 #no-penetration (wall)

def getAFBC_u(x,flag):
    if flag in exterior_boundary or flag == boundaryTags['right_porous']:
        return None #Dirichlet or outflow

    
def getAFBC_v(x,flag):
    if flag in exterior_boundary or flag == boundaryTags['right_porous']:
        return None #Dirichlet or outflow
   
advectiveFluxBoundaryConditions =  {0:getAFBC_p,
                                    1:getAFBC_u,
                                    2:getAFBC_v}
def getDFBC_u(x,flag):
    if flag in [boundaryTags['right_fluid']]:#outflow
        return lambda x,t: 0.0
    else: #Dirichlet
        return None

def getDFBC_v(x,flag):
    if flag in [boundaryTags['right_fluid']]:#outflow
        return lambda x,t: 0.0
    else: #Dirichlet
        return None

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u},
                                   2:{2:getDFBC_v}}

class Steady_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return -(L[1]-x[1])*coefficients.rho*coefficients.g[1]

class Steady_u:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

class Steady_v:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:Steady_p(),
                     1:Steady_u(),
                     2:Steady_v()}
## @}
