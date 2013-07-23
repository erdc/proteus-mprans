from proteus.default_so import *
import wavetank

if wavetank.useOnlyVF:
    pnList = [("twp_navier_stokes_p", #0
               "twp_navier_stokes_n"),
              ("vof_p", #1              
               "vof_n")]
else:
    pnList = [("twp_navier_stokes_p" , #0
               "twp_navier_stokes_n"),
              ("vof_p" , #1
               "vof_n"),
              ("ls_p" , #2
               "ls_n"),
              ("redist_p" ,#3 
               "redist_n"),
              ("ls_consrv_p" ,#4 
               "ls_consrv_n")]

if wavetank.movingDomain:
    pnList.append(("moveMesh_p","moveMesh_n"))

if wavetank.useRANS > 0:
    pnList.append(("kappa_p",
                   "kappa_n"))
    pnList.append(("dissipation_p",
                   "dissipation_n"))
name = "wavetank"

#systemStepControllerType = ISO_fixed_MinAdaptiveModelStep
systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,wavetank.dt_init]+[wavetank.dt_init+ i*wavetank.dt_out for i in range(1,wavetank.nDTout+1)]

