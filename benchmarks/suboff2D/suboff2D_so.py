from proteus.default_so import *
from suboff2D import *

pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
          ("ls_p",                "ls_n"),
          ("redist_p",            "redist_n")]
          #("kappa_p", "kappa_n"),
          #("dissipation_p", "dissipation_n")]

    
name = "suboff2D_p"

systemStepControllerType = Sequential_MinAdaptiveModelStep
    
needEBQ_GLOBAL = False
needEBQ = False


