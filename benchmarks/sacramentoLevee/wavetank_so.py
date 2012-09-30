from proteus.default_so import *
import wavetank

pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
          ("ls_p",                "ls_n"),
          ("vof_p",               "vof_n"),
          ("redist_p",            "redist_n"),
          ("ls_consrv_p",         "ls_consrv_n")]

    
name = "wavetank" 

systemStepControllerType = Sequential_MinAdaptiveModelStep#Sequential_FixedStep
    
needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,wavetank.dt_init]+[i*wavetank.dt_fixed for i in range(1,wavetank.nDTout+1)] 
