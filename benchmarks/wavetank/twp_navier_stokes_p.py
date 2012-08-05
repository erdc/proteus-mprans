from proteus import *
from proteus.default_p import *
from wavetank import *
from proteus.mprans import RANS2P

LevelModelType = RANS2P.LevelModel
coefficients = RANS2P.Coefficients(epsFact=epsFact_viscosity,
                                             sigma=0.0,
                                             rho_0 = rho_0,
                                             nu_0 = nu_0,
                                             rho_1 = rho_1,
                                             nu_1 = nu_1,
                                             g=g,
                                             nd=nd,
                                             LS_model=1,
                                             epsFact_density=epsFact_density,
                                             stokes=False,
                                             useRBLES=useRBLES,
		                             useMetrics=useMetrics)
					     

def getDBC_p(x,flag):
    if x[2] > L[2] - 1.0e-8:#top atmospheric
        return lambda x,t: rho_1*x[2]*g[2]
    elif x[0] > L[0] - 1.0e-8: #right hydrostatic w.r.t outflowHeight
        return outflowPressure

def getDBC_u(x,flag):
    if x[0] < 1.0e-8:
        return inflowVelocity_u
    elif x[2] > L[2] - 1.0e-8:
        return lambda x,t: windspeed_u

def getDBC_v(x,flag):
    if x[0] < 1.0e-8:
        return inflowVelocity_v
    elif x[2] > L[2] - 1.0e-8:
        return lambda x,t: windspeed_v

def getDBC_w(x,flag):
    if x[0] < 1.0e-8:
        return inflowVelocity_w
    
dirichletConditions = {0:getDBC_p,
                       1:getDBC_u,
                       2:getDBC_v,
                       3:getDBC_w}

def getAFBC_p(x,flag):
    if x[0] < 1.0e-8:
        return inflowFlux
    elif x[0] > L[0] - 1.0e-8: #right open
        return None
    elif x[2] > L[2] - 1.0e-8: #top open
        return None
    else:  #wall
        return lambda x,t: 0.0

def getAFBC_u(x,flag):
    if x[0] < 1.0e-8:#left Dirichlet
        return None
    elif x[0] > L[0] - 1.0e-8: #right open
        return None
    elif x[2] > L[2] - 1.0e-8: #top Dirichlet on x-component
        return None
    else:
        return lambda x,t: 0.0

def getAFBC_v(x,flag):
    if x[0] < 1.0e-8:#left Dirichlet
        return None
    elif x[0] > L[0] - 1.0e-8: #right outflow
        return None
    elif x[2] > L[2] - 1.0e-8: #top Dirichlet on y-component
        return None
    else:
        return lambda x,t: 0.0

def getAFBC_w(x,flag):
    if x[0] < 1.0e-8:#left Dirichlet
        return None
    elif x[0] > L[0] - 1.0e-8: #right outflow
        return None
    elif x[2] > L[2] - 1.0e-8: #top open
        return None
    else:
        return lambda x,t: 0.0

def getDFBC_u(x,flag):
    if x[0] < 1.0e-8:#left Dirichlet
        return None
    elif x[0] > L[0] - 1.0e-8: #right no diffusive flux (open)
        return lambda x,t: 0.0
    elif x[2] > L[2] - 1.0e-8: #top Dirichlet on x-component
        return None
    else:  #no flux everywhere else
        return lambda x,t: 0.0

def getDFBC_v(x,flag):
    if x[0] < 1.0e-8:#left Dirichlet
        return None
    elif x[0] > L[0] - 1.0e-8: #right outflow
        return lambda x,t: 0.0
    elif x[2] > L[2] - 1.0e-8: #top Dirichlet on y-component
        return None
    else:  #no flux everywhere else
        return lambda x,t: 0.0

def getDFBC_w(x,flag):
    if x[0] < 1.0e-8:#left Dirichlet
        return None
    elif x[2] > L[2] - 1.0e-8: #top outflow
        return lambda x,t: 0.0
    elif x[0] > L[0] - 1.0e-8: #right outflow
        return lambda x,t: 0.0
    else: #no diffusive flux everywhere else
        return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_p,
                                    1:getAFBC_u,
                                    2:getAFBC_v,
                                    3:getAFBC_w}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u},
                                   2:{2:getDFBC_v},
                                   3:{3:getDFBC_w}}

class PerturbedSurface_p:
    def __init__(self,waterLevel):
        self.waterLevel=waterLevel
    def uOfXT(self,x,t):
        return smoothedHydrostaticPressure(self.waterLevel,x[2])

class MeanFlow:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return inflowVelocityMean[0]

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:PerturbedSurface_p(inflowHeightMean),
                     1:MeanFlow(),
                     2:AtRest(),
                     3:AtRest()}
