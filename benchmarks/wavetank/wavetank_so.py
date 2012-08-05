from proteus.default_so import *
import wavetank

pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
          ("ls_p",                "ls_n"),
          ("vof_p",               "vof_n"),
          ("redist_p",            "redist_n"),
          ("ls_consrv_p",         "ls_consrv_n")]

    
name = "wavetank" 

systemStepControllerType = Sequential_FixedStep
    
needEBQ_GLOBAL = False
needEBQ = False

tnList = [i*wavetank.dt_fixed for i in range(0,wavetank.nDTout+1)] 
