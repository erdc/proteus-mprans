from proteus import *
from proteus.default_p import *
from wigley import *
from proteus.mprans import NCLSV2

LevelModelType = NCLSV2.OneLevelNCLSV2

coefficients = NCLevelSetCoefficients(V_model=0,RD_model=3,ME_model=1,checkMass=False)
#coefficients = NCLevelSetCoefficients(V_model=0,ME_model=1)

def getDBC_ls(x,flag):
    pass

dirichletConditions = {0:getDBC_ls}

class Flat_phi:
    def __init__(self,waterLevel):
        self.waterLevel=waterLevel
    def uOfXT(self,x,t):
        signedDistance = x[2] - self.waterLevel
        return signedDistance

initialConditions  = {0:Flat_phi(waterLevel)}
    
fluxBoundaryConditions = {0:'outFlow'}
    
advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
