from proteus.default_so import *
from suboff import *

pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
          ("ls_p",                "ls_n"),
          ("redist_p",            "redist_n"),
          ("kappa_p", "kappa_n"),
          ("epsilon_p", "epsilon_n")]

    
name = "suboff_p"

systemStepControllerType = Sequential_MinAdaptiveModelStep
    
needEBQ_GLOBAL = False
needEBQ = False


