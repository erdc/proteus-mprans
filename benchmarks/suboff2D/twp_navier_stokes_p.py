from proteus import *
from proteus.default_p import *
from suboff2D import *
from proteus.mprans import RANS2P

LevelModelType = RANS2P.LevelModel

turbulenceClosureModel = 3 #K-Epsilon
if dissipation_model_flag >= 2:
    turbulenceClosureModel = 4 #K-Omega

coefficients = RANS2P.Coefficients(epsFact=epsFact_viscosity,
                                   sigma=0.0,
                                   rho_0 = rho_0,
                                   nu_0 = nu_0,
                                   rho_1 = rho_0,
                                   nu_1 = nu_0,
                                   g=g,
                                   nd=nd,
                                   LS_model=1,
                                   Closure_0_model=None,#3,
                                   Closure_1_model=None,#4,
                                   turbulenceClosureModel=1,
                                   epsFact_density=epsFact_density,
                                   stokes=False,
                                   useRBLES=useRBLES,
                                   useMetrics=useMetrics,
                                   useVF=1.0,
                                   eb_adjoint_sigma=1.0,
                                   eb_penalty_constant=100.,
                                   forceStrongDirichlet=True)

getDBC_p = twpflowPressure
getDBC_u = twpflowVelocity_u
getDBC_v = twpflowVelocity_v

dirichletConditions = {0:getDBC_p,
                       1:getDBC_u,
                       2:getDBC_v}

def getAFBC_p(x,flag):
    if flag == boundaryTags['bottom'] or flag == boundaryTags['top'] or flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0

def getAFBC_u(x,flag):
    if flag == boundaryTags['bottom'] or flag == boundaryTags['top']:
        return lambda x,t: 0.0

def getAFBC_v(x,flag):
    if flag == boundaryTags['bottom'] or flag == boundaryTags['top']:
        return lambda x,t: 0.0

def getDFBC_u(x,flag):
    if not (flag == boundaryTags['left'] or flag == boundaryTags['obstacle']):
        return lambda x,t: 0.0
    
def getDFBC_v(x,flag):
    if not (flag == boundaryTags['left'] or flag == boundaryTags['obstacle']):
        return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_p,
                                    1:getAFBC_u,
                                    2:getAFBC_v}

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
        return inflowVelocity[0]

class Steady_v:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return inflowVelocity[1]
    

initialConditions = {0:Steady_p(),
                     1:Steady_u(),
                     2:Steady_v()}

## @}
