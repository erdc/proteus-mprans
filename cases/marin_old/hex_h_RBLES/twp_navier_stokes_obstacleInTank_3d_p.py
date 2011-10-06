from proteus import *
from proteus.default_p import *
from obstacleInTank3d import *
from proteus.mprans import RANS2P2
from proteus.ctransportCoefficients import smoothedHeaviside_integral

LevelModelType = RANS2P2.LevelModel

coefficients = RANS2P2.Coefficients(epsFact=epsFact_viscosity,
                                   sigma=0.0,
                                   rho_0 = rho_0,
                                   nu_0 = nu_0,
                                   rho_1 = rho_1,
                                   nu_1 = nu_1,
                                   g=g,
                                   nd=nd,
                                   LS_model=3,
                                   epsFact_density=epsFact_density,
                                   stokes=useStokes,
                                   useRBLES=useRBLES,
		                   useMetrics=useMetrics)
				   
				   

def getDBC_p_obstacleInTank(x,flag):
    return None
def getDBC_u_obstacleInTank(x,flag):
    return None
def getDBC_v_obstacleInTank(x,flag):
    return None
def getDBC_w_obstacleInTank(x,flag):
    return None

dirichletConditions = {0:getDBC_p_obstacleInTank,
                       1:getDBC_u_obstacleInTank,
                       2:getDBC_v_obstacleInTank,
                       3:getDBC_w_obstacleInTank}

eps=1.0e-4

def getAFBC_p_obstacleInTank(x,flag):
    if x[0] in [0.0, 3.22] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]:
        return lambda x,t: 0.0
    if x[0]>=2.3955-eps and x[0]<=2.5565+eps and x[1]>=0.2985-eps and x[1]<=0.7015+eps and x[2]>=0.0-eps and x[2]<=0.161+eps:
        return lambda x,t: 0.0
def getAFBC_u_obstacleInTank(x,flag):
    if x[0] in [0.0, 3.22] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]:
        return lambda x,t: 0.0
    if x[0]>=2.3955-eps and x[0]<=2.5565+eps and x[1]>=0.2985-eps and x[1]<=0.7015+eps and x[2]>=0.0-eps and x[2]<=0.161+eps:
        return lambda x,t: 0.0
def getAFBC_v_obstacleInTank(x,flag):
    if x[0] in [0.0, 3.22] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]:
        return lambda x,t: 0.0
    if x[0]>=2.3955-eps and x[0]<=2.5565+eps and x[1]>=0.2985-eps and x[1]<=0.7015+eps and x[2]>=0.0-eps and x[2]<=0.161+eps:
        return lambda x,t: 0.0
def getAFBC_w_obstacleInTank(x,flag):
    if x[0] in [0.0, 3.22] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]:
        return lambda x,t: 0.0
    if x[0]>=2.3955-eps and x[0]<=2.5565+eps and x[1]>=0.2985-eps and x[1]<=0.7015+eps and x[2]>=0.0-eps and x[2]<=0.161+eps:
        return lambda x,t: 0.0
#def getDFBC_u_obstacleInTank(x,flag):
#    if x[0] in [0.0, 3.22] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]:
#        return lambda x,t: 0.0
#def getDFBC_v_obstacleInTank(x,flag):
#    if x[0] in [0.0, 3.22] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]:
#        return lambda x,t: 0.0
#def getDFBC_w_obstacleInTank(x,flag):
#    if x[0] in [0.0, 3.22] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]:
#        return lambda x,t: 0.0

def getDFBC_u_obstacleInTank(x,flag):
    return lambda x,t: 0.0
def getDFBC_v_obstacleInTank(x,flag):
    return lambda x,t: 0.0
def getDFBC_w_obstacleInTank(x,flag):
    return lambda x,t: 0.0

fluxBoundaryConditions = {0:'mixedFlow',
                          1:'mixedFlow',
                          2:'mixedFlow',
                          3:'mixedFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_obstacleInTank,
                                    1:getAFBC_u_obstacleInTank,
                                    2:getAFBC_v_obstacleInTank,
                                    3:getAFBC_w_obstacleInTank}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u_obstacleInTank},
                                   2:{2:getDFBC_v_obstacleInTank},
                                   3:{3:getDFBC_w_obstacleInTank}}

class Shock_p:
    def uOfXT(self,x,t):
        if shockSignedDistance(x) < 0:
            return -g[2]*(rho_0*(height - x[2])
                          -(rho_0-rho_1)*smoothedHeaviside_integral(epsFact_density*0.1,height-waterColumn_z)
                          +(rho_0-rho_1)*smoothedHeaviside_integral(epsFact_density*0.1,x[2]-waterColumn_z))
        else:
            return -rho_1*g[2]

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:Shock_p(),
                     1:AtRest(),
                     2:AtRest(),
                     3:AtRest()}
