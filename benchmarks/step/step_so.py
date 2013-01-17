from proteus.default_so import *
from step import *

pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
          ("ls_p",                "ls_n"),
          ("redist_p",            "redist_n")]
#,
#          ("ls_p",                "ls_n"),
#          ("vof_p",               "vof_n"),
#          ("redist_p",            "redist_n"),
#          ("ls_consrv_p",         "ls_consrv_n")]

    
name = "step_p"

systemStepControllerType = Sequential_FixedStep
    
needEBQ_GLOBAL = False
needEBQ = False


