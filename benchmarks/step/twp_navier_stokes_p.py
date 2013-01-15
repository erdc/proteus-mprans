from proteus import *
from proteus.default_p import *
from step import *
from proteus.mprans import RANS2P

LevelModelType = RANS2P.LevelModel
coefficients = RANS2P.Coefficients(epsFact=epsFact_viscosity,
                                             sigma=0.0,
                                             rho_0 = rho_0,
                                             nu_0 = nu_0,
                                             rho_1 = rho_0,
                                             nu_1 = nu_0,
                                             g=g,
                                             nd=nd,
                                             LS_model=1,
                                             epsFact_density=epsFact_density,
                                             stokes=False,
                                             useRBLES=useRBLES,
		                             useMetrics=useMetrics)

bcsTimeDependent = True

periodic = False

getDBC_p = twpflowPressure
getDBC_u = twpflowVelocity_u
getDBC_v = twpflowVelocity_v
getDBC_w = twpflowVelocity_w

dirichletConditions = {0:getDBC_p,
                       1:getDBC_u,
                       2:getDBC_v,
                       3:getDBC_w}

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow',
                          3:'outFlow'}

def getAFBC_p(x,flag):
    if flag == 0:
        return lambda x,t: 0.0
    if flag == boundaryTags['upstream']:
        return lambda x,t: -inflow*velRamp(t)
    if (flag in [boundaryTags['front'],boundaryTags['back'],boundaryTags['top']]+bottom):
        return lambda x,t: 0.0

def getAFBC_u(x,flag):
    if flag == 0:
        return lambda x,t: 0.0
    if (flag in [boundaryTags['front'],boundaryTags['back']]):
        return lambda x,t: 0.0
    
def getAFBC_v(x,flag):
    if flag == 0:
        return lambda x,t: 0.0
    if (flag in [boundaryTags['front'],boundaryTags['back']]):
        return lambda x,t: 0.0
    
def getAFBC_w(x,flag):
    if flag == 0:
        return lambda x,t: 0.0
    if (flag in [boundaryTags['front'],boundaryTags['back']]):
        return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_p,
                                    1:getAFBC_u,
                                    2:getAFBC_v,
                                    3:getAFBC_w}
def getDFBC_u(x,flag):
    if flag == 0:
        return lambda x,t: 0.0
    if not (flag == boundaryTags['upstream'] or
            flag == boundaryTags['top'] or
            flag in bottom):
        return lambda x,t: 0.0

def getDFBC_v(x,flag):
    if flag == 0:
        return lambda x,t: 0.0
    if not (flag == boundaryTags['top'] or
            flag in bottom):
        return lambda x,t: 0.0

def getDFBC_w(x,flag):
    if flag == 0:
        return lambda x,t: 0.0
    if not (flag == boundaryTags['top'] or
            flag in bottom):
        return lambda x,t: 0.0

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u},
                                   2:{2:getDFBC_v},
                                   3:{3:getDFBC_w}}

class Steady_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return -(downstream_height-x[1])*coefficients.rho*coefficients.g[1]

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

class Steady_w:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:Steady_p(),
                     1:Steady_u(),
                     2:Steady_v(),
                     3:Steady_w()}

## @}
