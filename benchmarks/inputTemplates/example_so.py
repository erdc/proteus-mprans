"""
The split operator module for air/water flow around a moving rigid cylinder
"""
from proteus.default_so import *
import example

if example.useOnlyVF:
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

if example.movingDomain:
    pnList.append(("moveMesh_p","moveMesh_n"))

name = "example"

systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,example.dt_init]+[example.dt_init+ i*example.dt_out for i in range(1,example.nDTout+1)]

