from proteus import *
from proteus.default_p import *
from math import *
from beach_erosion_board_waves_3d import *
from proteus import RDLSV2
if useRDLS:
    LevelModelType = RDLSV2.OneLevelRDLSV2
"""
The redistancing equation in the sloshbox test problem.
"""
##

##\ingroup test
#\brief The redistancing equation in the sloshbox test problem.
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
if rd_freezeLS:
    if LevelModelType == RDLSV2.OneLevelRDLSV2:
        weakDirichletConditions = {0:RDLSV2.setZeroLSweakDirichletBCs}
    else:
        weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs}

initialConditions  = None

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
