from proteus import *
from proteus.default_p import *
from wigley import *
from proteus.mprans import NCLS

LevelModelType = NCLS.LevelModel

coefficients = NCLS.Coefficients(V_model=0,RD_model=3,ME_model=1,checkMass=checkMass,useMetrics=useMetrics,
				  sc_uref=ls_sc_uref,
				  sc_beta=ls_sc_beta)
#coefficients = NCLevelSetCoefficients(V_model=0,ME_model=1)

class Flat_phi:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return ls_wave(x,t)

analyticalSolutions = None

def getDBC_ls(x,flag):
    if flag == boundaryTags['left']:
        return ls_wave
	    
dirichletConditions = {0:getDBC_ls}

initialConditions  = {0:Flat_phi()}
    
fluxBoundaryConditions = {0:'mixedFlow'}

def getAFBC_ls(x,flag):
    return lambda x,t: 0.0
    
advectiveFluxBoundaryConditions =  {0:getAFBC_ls}

diffusiveFluxBoundaryConditions = {0:{}}
