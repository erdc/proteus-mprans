from proteus import *
from proteus.default_p import *
from obstacleInTank3d import *
from proteus import NCLS
from proteus import NCLSV2

LevelModelType = NCLS.OneLevelNCLS
LevelModelType = NCLSV2.OneLevelNCLSV2

if applyCorrection:
    coefficients = NCLevelSetCoefficients(V_model=0,RD_model=3,ME_model=1,checkMass=False)
elif applyRedistancing:
    coefficients = NCLevelSetCoefficients(V_model=0,RD_model=2,ME_model=1,checkMass=False)
else:
    coefficients = NCLevelSetCoefficients(V_model=0,ME_model=1,checkMass=False)

def getDBC_ls(x,flag):
    pass

dirichletConditions = {0:getDBC_ls}


class Shock_phi:
    def uOfXT(self,x,t):
        return shockSignedDistance(x)

initialConditions  = {0:Shock_phi()}
    
fluxBoundaryConditions = {0:'outFlow'}
    
advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
