from proteus import *
from proteus.default_p import *
from math import *
from obstacleInTank3d import *
from proteus import RDLS
from proteus import RDLSV2
"""
The redistancing equation in the obstacleInTank test problem.
"""
##
LevelModelType = RDLS.OneLevelRDLS
LevelModelType = RDLSV2.OneLevelRDLSV2

##\ingroup test
#\brief The redistancing equation in the obstacleInTank test problem.
#
if applyCorrection:
    coefficients = RedistanceLevelSet(applyRedistancing=applyRedistancing,
                                      epsFact=epsFact_redistance,
                                      nModelId=1,
                                      rdModelId=3)
else:
    coefficients = RedistanceLevelSet(applyRedistancing=applyRedistancing,
                                      epsFact=epsFact_redistance,
                                      nModelId=1,
                                      rdModelId=2)

#now define the Dirichlet boundary conditions

def getDBC_rd(x,flag):
    pass
    
dirichletConditions = {0:getDBC_rd}

if freezeLevelSet:
    if LevelModelType == RDLS.OneLevelRDLS:
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
