from proteus.default_so import *
import dambreak

pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
          ("ls_p",                "ls_n"),
          ("vof_p",               "vof_n"),
          ("redist_p",            "redist_n"),
          ("ls_consrv_p",         "ls_consrv_n")]

    
name = "dambreak_p" 

systemStepControllerType = Sequential_MinAdaptiveModelStep#Sequential_FixedStep
    
needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,dambreak.dt_init]+[i*dambreak.dt_fixed for i in range(1,dambreak.nDTout+1)] 

