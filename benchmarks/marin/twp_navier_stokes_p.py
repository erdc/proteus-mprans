from proteus import *
from proteus.default_p import *
from marin import *
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
    return None

def getDBC_u(x,flag):
    return None

def getDBC_v(x,flag):
    return None

def getDBC_w(x,flag):
    return None
    
dirichletConditions = {0:getDBC_p,
                       1:getDBC_u,
                       2:getDBC_v,
                       3:getDBC_w}

def getAFBC_p(x,flag):
    return lambda x,t: 0.0

def getAFBC_u_old(x,flag):
    if (flag == bt['left']     or
        flag == bt['right']    or
	flag == bt['box_left'] or
        flag == bt['box_right']):
       return lambda x,t: 0.0

def getAFBC_v_old(x,flag):
    if (flag == bt['front']     or
        flag == bt['back']      or
	flag == bt['box_front'] or
        flag == bt['box_back'] ):
       return lambda x,t: 0.0

def getAFBC_w_old(x,flag):
    if (flag == bt['bottom']     or
        flag == bt['top']        or
        flag == bt['box_top']):
       return lambda x,t: 0.0

eps=1.0e-4
     
def getAFBC_u(x,flag):
    if x[0] in [0.0, 3.22] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]:
        return lambda x,t: 0.0
    if x[0]>=2.3955-eps and x[0]<=2.5565+eps and x[1]>=0.2985-eps and x[1]<=0.7015+eps and x[2]>=0.0-eps and x[2]<=0.161+eps:
        return lambda x,t: 0.0
	
def getAFBC_v(x,flag):
    if x[0] in [0.0, 3.22] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]:
        return lambda x,t: 0.0
    if x[0]>=2.3955-eps and x[0]<=2.5565+eps and x[1]>=0.2985-eps and x[1]<=0.7015+eps and x[2]>=0.0-eps and x[2]<=0.161+eps:
        return lambda x,t: 0.0
	
def getAFBC_w(x,flag):
    if x[0] in [0.0, 3.22] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]:
        return lambda x,t: 0.0
    if x[0]>=2.3955-eps and x[0]<=2.5565+eps and x[1]>=0.2985-eps and x[1]<=0.7015+eps and x[2]>=0.0-eps and x[2]<=0.161+eps:
        return lambda x,t: 0.0
	

def getDFBC_u(x,flag):
    return lambda x,t: 0.0
    
def getDFBC_v(x,flag):
    return lambda x,t: 0.0

def getDFBC_w(x,flag):
    return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_p,
                                    1:getAFBC_u,
                                    2:getAFBC_v,
                                    3:getAFBC_w}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u},
                                   2:{2:getDFBC_v},
                                   3:{3:getDFBC_w}}

fluxBoundaryConditions = {0:'mixedFlow',
                          1:'mixedFlow',
                          2:'mixedFlow',
                          3:'mixedFlow'}

class PerturbedSurface_p:
    def __init__(self,waterLevel):
        self.waterLevel=waterLevel
    def uOfXT(self,x,t):
        if signedDistance(x) < 0:
            return -(L[2] - self.waterLevel)*rho_1*g[2] - (self.waterLevel - x[2])*rho_0*g[2]
        else:
            return -(L[2] - self.waterLevel)*rho_1*g[2]

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:PerturbedSurface_p(waterLine_z),
                     1:AtRest(),
                     2:AtRest(),
                     3:AtRest()}
