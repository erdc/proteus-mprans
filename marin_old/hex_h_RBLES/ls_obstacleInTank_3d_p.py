from proteus import *
from proteus.default_p import *
from obstacleInTank3d import *
from proteus.mprans import NCLS

LevelModelType = NCLS.LevelModel

coefficients = NCLS.Coefficients(V_model=0,RD_model=3,ME_model=1,checkMass=False,useMetrics=useMetrics)


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
