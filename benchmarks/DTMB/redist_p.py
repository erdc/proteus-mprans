from proteus import *
from proteus.default_p import *
from math import *
from DTMB import *
from proteus.mprans import RDLS
"""
The redistancing equation in the sloshbox test problem.
"""

LevelModelType = RDLS.LevelModel

coefficients = RDLS.Coefficients(applyRedistancing=True,
                                 epsFact=epsFact_redistance,
                                 nModelId=1,
                                 rdModelId=3,
		                 useMetrics=useMetrics)

def getDBC_rd(x,flag):
    pass
    
dirichletConditions     = {0:getDBC_rd}
weakDirichletConditions = {0:RDLS.setZeroLSweakDirichletBCsSimple}

advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{}}

class PHI_IC:       
    def uOfXT(self,x,t):
        return wavePhi(x,t)

initialConditions  = {0:PHI_IC()}
