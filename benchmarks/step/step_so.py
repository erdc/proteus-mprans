from proteus.default_so import *
from step import *

pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
          ("ls_p",                "ls_n"),
          ("redist_p",            "redist_n"),
          ("kappa_p", "kappa_n"),
          ("dissipation_p", "dissipation_n")]

    
name = "step_p"

systemStepControllerType = Sequential_MinAdaptiveModelStep
    
needEBQ_GLOBAL = False
needEBQ = False


