from proteus import *
from proteus.default_p import *
from beach_erosion_board_waves_3d import *
from proteus.mprans import VOFV2
from proteus import VolumeAveragedVOF

if useVOF:
    if useSpongeLayer:
        LevelModelType = VOFV2.OneLevelVOFV2
    else:
        LevelModelType = VolumeAveragedVOF.OneLevelVolumeAveragedVOF

if useSpongeLayer:
    coefficients = VolumeAveragedVOFCoefficients(LS_model=1,V_model=0,RD_model=3,ME_model=2,epsFact=epsFact_vof,
                                                 setParamsFunc=spongeLayerFunc,checkMass=checkMass)
else:
    coefficients = VOFCoefficients(LS_model=1,V_model=0,RD_model=3,ME_model=2,epsFact=epsFact_vof,checkMass=checkMass)

dirichletConditions = {0:getDBC_vof}

initialConditions  = {0:Flat_H()}

fluxBoundaryConditions = {0:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_vof}

diffusiveFluxBoundaryConditions = {0:{}}
