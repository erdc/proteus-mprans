from proteus.default_so import *
import sloshbox3d
from sloshbox3d import *

if sloshbox3d.applyCorrection:
    pnList = [("twp_navier_stokes_sloshbox_3d_p" , 
               "twp_navier_stokes_sloshbox_3d_n"),
              ("ls_sloshbox_3d_p" , "ls_sloshbox_3d_n"),
              ("vof_sloshbox_3d_p" , "vof_sloshbox_3d_n"),
              ("redist_sloshbox_3d_p" , "redist_sloshbox_3d_n"),
              ("ls_consrv_sloshbox_3d_p" , 
               "ls_consrv_sloshbox_3d_n")]
elif sloshbox3d.applyRedistancing:
    pnList = [("twp_navier_stokes_sloshbox_3d_p" , 
               "twp_navier_stokes_sloshbox_3d_n"),
              ("ls_sloshbox_3d_p" , "ls_sloshbox_3d_n"),
              ("redist_sloshbox_3d_p" , "redist_sloshbox_3d_n")]
else:
    pnList = [("twp_navier_stokes_sloshbox_3d_p" , 
               "twp_navier_stokes_sloshbox_3d_n"),
              ("ls_sloshbox_3d_p" , "ls_sloshbox_3d_n")]
    
name = "twp_navier_stokes_sloshbox_3d"

systemStepControllerType = Sequential_FixedStep
    
needEBQ_GLOBAL = False#True
needEBQ = False#True
useOneArchive = True
name=soname


tnList = [i*sloshbox3d.dt_fixed for i in range(0,sloshbox3d.nDTout+1)] 

#tnList = [0.0,sloshbox3d.dt_init,1.0,2.0,10.0,20.0]
