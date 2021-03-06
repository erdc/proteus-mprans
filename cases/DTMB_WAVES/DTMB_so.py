"""
The split operator module for air/water flow around a moving rigid cylinder
"""
from proteus.default_so import *
import DTMB

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

name = "DTMB"

systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False
useOneArchive = False
archiveFlag = ArchiveFlags.EVERY_USER_STEP
tnList = [0.0,DTMB.dt_init]+[DTMB.dt_init*i for i in range(2,DTMB.nDTout,1)]

