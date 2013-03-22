from proteus import *
from proteus.default_p import *
from flume import *
from proteus.mprans import RANS2P

LevelModelType = RANS2P.LevelModel
if useOnlyVF:
    LS_model = None
else:
    LS_model = 2
coefficients = RANS2P.Coefficients(epsFact=epsFact_viscosity,
                                   sigma=0.0,
                                   rho_0 = rho_0,
                                   nu_0 = nu_0,
                                   rho_1 = rho_1,
                                   nu_1 = nu_1,
                                   g=g,
                                   nd=nd,
                                   VF_model=1,
                                   LS_model=LS_model,
                                   epsFact_density=epsFact_density,
                                   stokes=False,
                                   useVF=useVF,
				   useRBLES=useRBLES,
				   useMetrics=useMetrics,
                                   eb_adjoint_sigma=1.0,
                                   forceStrongDirichlet=ns_forceStrongDirichlet,
                                   turbulenceClosureModel=ns_closure)

def getDBC_p(x,flag):
    if flag == boundaryTags['right']:
        return outflowPressure

def getDBC_u(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: inflowVelocity[0]
    if flag == boundaryTags['right']:
        return lambda x,t: outflowVelocity[0]
    if flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0

def getDBC_v(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: inflowVelocity[1]
    if flag == boundaryTags['right']:
        return lambda x,t: outflowVelocity[1]
    if flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0

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

class Hydrostatic_p:
    def uOfXT(self,x,t):
        return inflowPressure(x,t)

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

initialConditions = {0:Hydrostatic_p(),
                     1:Steady_u(),
                     2:Steady_v()}
