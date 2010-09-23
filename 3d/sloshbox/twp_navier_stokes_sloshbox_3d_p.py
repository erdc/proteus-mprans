from proteus import *
from proteus.default_p import *
from sloshbox3d import *
from proteus import RANS2PV2

useOpt=True#False
if useOpt:
    LevelModelType = RANS2PV2.OneLevelRANS2PV2
coefficients = TwophaseNavierStokes_ST_LS_SO(epsFact=epsFact_viscosity,
                                             sigma=0.0,
                                             rho_0 = rho_0,
                                             nu_0 = nu_0,
                                             rho_1 = rho_1,
                                             nu_1 = nu_1,
                                             g=g,
                                             nd=nd,
                                             LS_model=1,
                                             epsFact_density=epsFact_density,
                                             stokes=useStokes)
