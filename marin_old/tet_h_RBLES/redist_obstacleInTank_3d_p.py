from proteus import *
from proteus.default_p import *
from math import *
from obstacleInTank3d import *
from proteus.mprans import RDLS
"""
The redistancing equation in the obstacleInTank test problem.
"""
##
LevelModelType = RDLS.LevelModel

##\ingroup test
#\brief The redistancing equation in the obstacleInTank test problem.
#
coefficients = RDLS.Coefficients(applyRedistancing=True,
                                 epsFact=epsFact_redistance,
                                 nModelId=1,
                                 rdModelId=3,
		                 useMetrics=useMetrics)


#now define the Dirichlet boundary conditions

def getDBC_rd(x,flag):
    pass
    
dirichletConditions = {0:getDBC_rd}

if freezeLevelSet:
    if LevelModelType == RDLS.LevelModel:
        weakDirichletConditions = {0:RDLS.setZeroLSweakDirichletBCs}
    else:
        weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs}

class Shock_phi:
    def uOfXT(self,x,t):
        return shockSignedDistance(x)

initialConditions  = {0:Shock_phi()}

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
