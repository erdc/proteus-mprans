from proteus.default_so import *
from windtunnel import *

pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
          ("ls_p",                "ls_n"),
          ("redist_p",            "redist_n")]
if useRANS > 0:
    pnList.append(("kappa_p",
                   "kappa_n"))
    pnList.append(("dissipation_p",
                   "dissipation_n"))

name = "windtunnel_p"

systemStepControllerType = Sequential_MinAdaptiveModelStep
    
needEBQ_GLOBAL = False
needEBQ = False


