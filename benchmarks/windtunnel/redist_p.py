from proteus import *
from proteus.default_p import *
from math import *
from windtunnel import *
from proteus.mprans import RDLS
"""
The redistancing equation in the sloshbox test problem.
"""

LevelModelType = RDLS.LevelModel

coefficients = RDLS.Coefficients(applyRedistancing=False,
                                 epsFact=epsFact_redistance,
                                 nModelId=1,
                                 rdModelId=2,
		                 useMetrics=useMetrics)

def getDBC_rd(x,flag):
    pass
    
dirichletConditions     = {0:getDBC_rd}
weakDirichletConditions = {0:RDLS.setZeroLSweakDirichletBCsSimple}

advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{}}

class PerturbedSurface_phi:       
    def uOfXT(self,x,t):
        return signedDistance(x)
    
initialConditions  = {0:PerturbedSurface_phi()}
