from proteus.default_so import *
import marin
from marin import *

pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
          ("ls_p",                "ls_n"),
          ("vof_p",               "vof_n"),
          ("redist_p",            "redist_n"),
          ("ls_consrv_p",         "ls_consrv_n")]

    
name = "marin_p"

systemStepControllerType = Sequential_FixedStep
    
needEBQ_GLOBAL = False
needEBQ = False

tnList = [i*marin.dt_fixed for i in range(0,marin.nDTout+1)] 

