from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from wigley import *
from proteus.mprans import VOF

LevelModelType = VOF.LevelModel
coefficients = VOF.Coefficients(LS_model=1,
                               V_model=0,
                               RD_model=3,
                               ME_model=2,
                               epsFact=epsFact_vof,
                               checkMass=checkMass,
			       useMetrics=useMetrics,
                               sc_uref=vof_sc_uref,sc_beta=vof_sc_beta)

class Flat_H:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return vof_wave(x,t)

def getDBC_vof(x,flag):
    if flag == boundaryTags['left']:
        return vof_wave

def getAFBC_vof(x,flag):
    return lambda x,t: 0.0

    if [ flag == boundaryTags['front']  or
         flag == boundaryTags['back']   or
         flag == boundaryTags['top']    or
         flag == boundaryTags['bottom'] or
         flag == boundaryTags['right'] ]:
        return lambda x,t: 0.0


initialConditions  = {0:Flat_H()}
dirichletConditions = {0:getDBC_vof}
advectiveFluxBoundaryConditions =  {0:getAFBC_vof}
diffusiveFluxBoundaryConditions = {0:{}}

analyticalSolutions = None
#fluxBoundaryConditions = {0:'noFlow'}
fluxBoundaryConditions = {0:'mixedFlow'}
