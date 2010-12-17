from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from obstacleInTank3d import *
from proteus.mprans import VOF

LevelModelType = VOF.LevelModel
coefficients = VOF.Coefficients(LS_model=1,V_model=0,RD_model=3,ME_model=2,epsFact=epsFact_vof,checkMass=False,useMetrics=useMetrics,
                                sc_uref=vof_sc_uref,sc_beta=vof_sc_beta)

class Shock_H:
    def uOfXT(self,x,t):
        return smoothedHeaviside(epsFact_consrv_heaviside*0.1,shockSignedDistance(x))

analyticalSolutions = None

def getDBC_vof(x,flag):
    pass

dirichletConditions = {0:getDBC_vof}

initialConditions  = {0:Shock_H()}

fluxBoundaryConditions = {0:'mixedFlow'}

def getAFBC_vof(x,flag):
    return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_vof}

diffusiveFluxBoundaryConditions = {0:{}}
