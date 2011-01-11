from proteus import *
from proteus.default_p import *
from wigley import *
from proteus.mprans import NCLS

LevelModelType = NCLS.LevelModel

coefficients = NCLS.Coefficients(V_model=0,RD_model=3,ME_model=1,checkMass=checkMass,useMetrics=useMetrics,
				  sc_uref=ls_sc_uref,
				  sc_beta=ls_sc_beta)
#coefficients = NCLevelSetCoefficients(V_model=0,ME_model=1)

def getDBC_ls(x,flag):
    def ls(x,t):
        return x[2]-waterLevel
    if flag == domain.boundaryTags['left']:
        return ls
    if flag == domain.boundaryTags['right']:
        return ls
    if openSides:
        if flag == domain.boundaryTags['front']:
            return ls
        if flag == domain.boundaryTags['back']:
            return ls

dirichletConditions = {0:getDBC_ls}

class Flat_phi:
    def __init__(self,waterLevel):
        self.waterLevel=waterLevel
    def uOfXT(self,x,t):
        signedDistance = x[2] - self.waterLevel
        return signedDistance

initialConditions  = {0:Flat_phi(waterLevel)}
    
fluxBoundaryConditions = {0:'mixedFlow'}
    
advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
