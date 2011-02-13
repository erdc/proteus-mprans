from proteus import *
from proteus.default_p import *
from math import *
from floatingCylinder import *
from proteus.mprans import RDLS

LevelModelType = RDLS.LevelModel
coefficients = RDLS.Coefficients(applyRedistancing=applyRedistancing,
                                 epsFact=epsFact_redistance,
                                 nModelId=1,
                                 rdModelId=3,
                                 useConstantH=useConstantH)
#now define the Dirichlet boundary conditions

def getDBC_rd(x,flag):
    pass
    
dirichletConditions = {0:getDBC_rd}

if freezeLevelSet:
    if LevelModelType == RDLS.LevelModel:
        weakDirichletConditions = {0:RDLS.setZeroLSweakDirichletBCs}
    else:
        weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs}

class PerturbedSurface_phi:
    def __init__(self,waterLevel,slopeAngle):
        self.waterLevel=waterLevel
        self.slopeAngle=slopeAngle
    def uOfXT(self,x,t):
        if useShock:
            return shockSignedDistance(x)
        else:
            surfaceNormal = [-sin(self.slopeAngle),0,cos(self.slopeAngle)]
            signedDistance = (x[0] - 0.5)*surfaceNormal[0]+(x[2] - self.waterLevel)*surfaceNormal[2]
            return signedDistance
    
initialConditions  = {0:PerturbedSurface_phi(waterLevel,0.0)}

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
