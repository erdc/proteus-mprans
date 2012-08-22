from proteus import *
from proteus.default_p import *
from beach_erosion_board_waves_3d import *
#from proteus import VolumeAveragedVOF
from proteus.mprans import VOF#VOFV2
#from proteus.mprans import VolumeAveragedVOF


if useSpongeLayer:
    if useVOF:
        LevelModelType = VOF.LevelModel#VOF.LevelModel
        coefficients = VOF.Coefficients(LS_model=1,V_model=0,RD_model=3,ME_model=2,epsFact=epsFact_vof,checkMass=checkMass)
else:
    coefficients = VOF.Coefficients(LS_model=1,V_model=0,RD_model=3,ME_model=2,epsFact=epsFact_vof,checkMass=checkMass)

dirichletConditions = {0:getDBC_vof}

initialConditions  = {0:Flat_H()}

fluxBoundaryConditions = {0:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_vof}

diffusiveFluxBoundaryConditions = {0:{}}
