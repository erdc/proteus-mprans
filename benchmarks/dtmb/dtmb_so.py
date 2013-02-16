"""
The split operator module for air/water flow around a moving rigid cylinder
"""
from proteus.default_so import *
import dtmb

pnList = [("twp_navier_stokes_p" , #0
           "twp_navier_stokes_n"),
          ("ls_p" , #1
           "ls_n"),
          ("vof_p" , #2
           "vof_n"),
          ("redist_p" ,#3 
           "redist_n"),
          ("ls_consrv_p" ,#4 
           "ls_consrv_n"),
          ("moveMesh_p",#5
           "moveMesh_n")]
pnList = [("twp_navier_stokes_p" , #0
           "twp_navier_stokes_n"),
          ("ls_p" , #1
           "ls_n"),
          ("vof_p" , #2
           "vof_n"),
          ("redist_p" ,#3 
           "redist_n"),
          ("ls_consrv_p" ,#4 
           "ls_consrv_n"),
          ("kappa_p",
           "kappa_n"),#5
          ("epsilon_p",
           "epsilon_n")] #6

name = "dtmb"

systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,dtmb.dt_init]+[dtmb.dt_init+ i*dtmb.dt_out for i in range(1,dtmb.nDTout+1)]

