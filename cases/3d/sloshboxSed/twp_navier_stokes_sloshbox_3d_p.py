from proteus import *
from proteus.default_p import *
from sloshbox3d import *
from proteus.mprans import RANS3P

useOpt=True#False
if useOpt:
    LevelModelType = RANS3P.LevelModel
coefficients = RANS3P.Coefficients(epsFact=epsFact_viscosity,
                                   sigma=0.0,
                                   rho_0 = rho_0,
                                   nu_0 = nu_0,
                                   rho_1 = rho_1,
                                   nu_1 = nu_1,
                                   g=g,
                                   nd=nd,
                                   SED_model_ID =  None,#no sed effect on fluid for now
                                   #SED_model_ID =  5,
                                   LS_model=1,
                                   epsFact_density=epsFact_density,
                                   stokes=useStokes)
