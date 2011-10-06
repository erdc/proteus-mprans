from proteus import *
from proteus.default_p import *
from beach_erosion_board_waves_3d import *
from math import *

nd = 3

#coefficients = LevelSetCurvatureCoefficients(epsFact=0.01,LSModel_index=2)
coefficients = LevelSetCurvatureCoefficients(epsFact=epsFact_curvature,LSModel_index=2,nd=nd)

fluxBoundaryConditions = {0:'outFlow'}

def getDBC_kappa(x,flag):
    pass

dirichletConditions = {0:getDBC_kappa}

def getAFBC_kappa(x,flag):
     pass
advectiveFluxBoundaryConditions =  {0:getAFBC_kappa}

diffusiveFluxBoundaryConditions = {0:{}}
## @}
