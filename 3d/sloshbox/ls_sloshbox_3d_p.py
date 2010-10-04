from proteus import *
from proteus.default_p import *
from sloshbox3d import *
from proteus.mprans import NCLS

LevelModelType = NCLS.LevelModel

if applyCorrection:
    coefficients = NCLS.Coefficients(V_model=0,RD_model=3,ME_model=1,epsFact=epsFact_consrv_heaviside,checkMass=checkMass)
elif applyRedistancing:
    coefficients = NCLS.Coefficients(V_model=0,RD_model=2,ME_model=1,epsFact=epsFact_consrv_heaviside,checkMass=checkMass)
else:
    coefficients = NCLS.Coefficients(V_model=0,ME_model=1,epsFact=epsFact_consrv_heaviside,checkMass=checkMass)

def getDBC_ls(x,flag):
    pass

dirichletConditions = {0:getDBC_ls}

class PerturbedSurface_phi:
    def __init__(self,waterLevel,slopeAngle):
        self.waterLevel=waterLevel
        self.slopeAngle=slopeAngle
    def uOfXT(self,x,t):
        if useShock:
            return shockSignedDistance(x)
        else:
            surfaceNormal = [-sin(self.slopeAngle),cos(self.slopeAngle)]
            signedDistance = (x[0] - 0.5)*surfaceNormal[0]+(x[2] - self.waterLevel)*surfaceNormal[1]
            return signedDistance
    
initialConditions  = {0:PerturbedSurface_phi(waterLevel,slopeAngle)}
    
fluxBoundaryConditions = {0:'outFlow'}
    
advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

