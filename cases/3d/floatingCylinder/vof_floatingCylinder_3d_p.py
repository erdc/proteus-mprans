from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from floatingCylinder import *
from proteus.mprans import VOF

LevelModelType = VOF.LevelModel
coefficients = VOF.Coefficients(LS_model=1,V_model=0,RD_model=3,ME_model=2,epsFact=epsFact_consrv_heaviside,checkMass=False)

class Flat_H:
    def __init__(self,waterLevel):
        self.waterLevel=waterLevel
    def uOfXT(self,x,t):
        return smoothedHeaviside(epsFact_consrv_heaviside*he,x[2]-waterLevel)

analyticalSolutions = None

def getDBC_vof(x,flag):
    if altBC:
        if flag in [domain.boundaryTags['downstream'],domain.boundaryTags['upstream'],domain.boundaryTags['top']]:
            return lambda x,t: smoothedHeaviside(epsFact_consrv_heaviside*he,x[2]-waterLevel)
    else:
        if flag in [domain.boundaryTags['upstream']]:
            return lambda x,t: smoothedHeaviside(epsFact_consrv_heaviside*he,x[2]-waterLevel)
        if flag in [domain.boundaryTags['downstream']]:
            return lambda x,t: smoothedHeaviside(epsFact_consrv_heaviside*he,x[2]-waterLevel)

dirichletConditions = {0:getDBC_vof}

initialConditions  = {0:Flat_H(waterLevel)}

fluxBoundaryConditions = {0:'mixedFlow'}

def getAFBC_vof(x,flag):
    if altBC:
        if flag not in [domain.boundaryTags['downstream'],domain.boundaryTags['upstream'],domain.boundaryTags['top']]:
            return lambda x,t: 0.0
    else:
        if flag not in [domain.boundaryTags['upstream'],domain.boundaryTags['downstream']]:
            return lambda x,t: 0.0
    if flag == 0:
        return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_vof}

diffusiveFluxBoundaryConditions = {0:{}}
