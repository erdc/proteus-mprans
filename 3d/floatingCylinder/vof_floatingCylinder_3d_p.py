from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from floatingCylinder import *
from proteus.mprans import VOF

LevelModelType = VOF.LevelModel
coefficients = VOF.Coefficients(LS_model=1,V_model=0,RD_model=3,ME_model=2,epsFact=epsFact_vof,checkMass=False)

class PerturbedSurface_H:
    def __init__(self,waterLevel,slopeAngle):
        self.waterLevel=waterLevel
        self.slopeAngle=slopeAngle
    def uOfXT(self,x,t):
        if useShock:
            return smoothedHeaviside(epsFact_consrv_heaviside*he,shockSignedDistance(x))
        else:
            surfaceNormal = [-sin(self.slopeAngle),0,cos(self.slopeAngle)]
            signedDistance = (x[0] - 0.5)*surfaceNormal[0]+(x[2] - self.waterLevel)*surfaceNormal[2]
            return smoothedHeaviside(epsFact_consrv_heaviside*he,signedDistance)

analyticalSolutions = None

def getDBC_vof(x,flag):
    if flag == domain.boundaryTags['top']:
        return lambda x,t: 1.0
    if flag == domain.boundaryTags['upstream']:
        return lambda x,t: smoothedHeaviside(epsFact_consrv_heaviside*he,x[2]-waterLevel)
    if flag == domain.boundaryTags['downstream']:
        return lambda x,t: smoothedHeaviside(epsFact_consrv_heaviside*he,x[2]-waterLevel)

dirichletConditions = {0:getDBC_vof}

initialConditions  = {0:PerturbedSurface_H(waterLevel,0.0)}

#fluxBoundaryConditions = {0:'noFlow'}
fluxBoundaryConditions = {0:'mixedFlow'}

def getAFBC_vof(x,flag):
    if flag == boundaryTags['bottom']:
        return lambda x,t: 0.0
    if flag == boundaryTags['front']:
        return lambda x,t: 0.0
    if flag == boundaryTags['back']:
        return lambda x,t: 0.0
    if flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0
    if flag == 0:
        return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_vof}

diffusiveFluxBoundaryConditions = {0:{}}
