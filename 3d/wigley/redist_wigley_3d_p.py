from proteus import *
from proteus.default_p import *
from math import *
from wigley import *
from proteus import RDLSV2

LevelModelType = RDLSV2.OneLevelRDLSV2
coefficients = RedistanceLevelSet(applyRedistancing=applyRedistancing,
                                  epsFact=epsFact_redistance,
                                  nModelId=1,
                                  rdModelId=3)

def getDBC_rd(x,flag):
    pass
    
dirichletConditions = {0:getDBC_rd}

if freezeLevelSet:
    if LevelModelType == RDLSV2.OneLevelRDLSV2:
        weakDirichletConditions = {0:RDLSV2.setZeroLSweakDirichletBCs}
    else:
        weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs}
#weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs2}

class Flat_phi:
    def __init__(self,waterLevel):
        self.waterLevel=waterLevel
    def uOfXT(self,x,t):
        signedDistance = x[2] - self.waterLevel
        return signedDistance

initialConditions  = {0:Flat_phi(waterLevel)}

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
