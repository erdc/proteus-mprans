from proteus.default_so import *
import vortex
from vortex import *

if applyRedistancing:
    if applyCorrection:
        pnList = [("ls_vortex_3d_p","ls_vortex_3d_n"),
                  ("redist_vortex_3d_p","redist_vortex_3d_n"), 
                  ("vof_vortex_3d_p","vof_vortex_3d_n"), 
                  ("ls_consrv_vortex_3d_p","ls_consrv_vortex_3d_n")]
    else:
        pnList = [("ls_vortex_3d_p","ls_vortex_3d_n"),
                  ("redist_vortex_3d_p","redist_vortex_3d_n"), 
                  ("vof_vortex_3d_p","vof_vortex_3d_n")]
else:
    pnList = [("ls_vortex_3d_p","ls_vortex_3d_n"),
              ("vof_vortex_3d_p","vof_vortex_3d_n")]

if "FLCBDF" in [timeIntegration_vof,timeIntegration_ls]:
    systemStepControllerType = Sequential_MinFLCBDFModelStep
else:
    systemStepControllerType = Sequential_MinAdaptiveModelStep

name=soname
if tryOpt:
    needEBQ_GLOBAL  = False
    needEBQ = False
else:
    needEBQ_GLOBAL  = True
    needEBQ = True

if tryOpt:
    nDTout = 21
    archiveFlag = ArchiveFlags.EVERY_USER_STEP
    
else:
    nDTout = 1
DT = T/float(nDTout)
tnList = [i*DT for i  in range(nDTout+1)]
#cek hard coded steps for article snapshots
#tnList = [0.0,4.0,8.0]
useOneArchive = True
