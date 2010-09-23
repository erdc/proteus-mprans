from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from obstacleInTank3d import *
from proteus import VOF
from proteus import VOFV2

LevelModelType = VOF.OneLevelVOF
LevelModelType = VOFV2.OneLevelVOFV2
coefficients = VOFCoefficients(LS_model=1,V_model=0,RD_model=3,ME_model=2,epsFact=epsFact_vof,checkMass=False)

class Shock_H:
    def uOfXT(self,x,t):
        return smoothedHeaviside(epsFact_consrv_heaviside*he,shockSignedDistance(x))

analyticalSolutions = None

def getDBC_vof(x,flag):
    if flag == boundaryTags['top']:
        return lambda x,t: 1.0

dirichletConditions = {0:getDBC_vof}

initialConditions  = {0:Shock_H()}

fluxBoundaryConditions = {0:'mixedFlow'}

def getAFBC_vof(x,flag):
    if flag != boundaryTags['top']:
        return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_vof}

diffusiveFluxBoundaryConditions = {0:{}}
