"""
The split operator module for air/water flow around a moving rigid cylinder
"""
from proteus.default_so import *
import floatingCylinder

pnList = [("twp_navier_stokes_floatingCylinder_3d_p" , #0
           "twp_navier_stokes_floatingCylinder_3d_n"),
          ("ls_floatingCylinder_3d_p" , #1
           "ls_floatingCylinder_3d_n"),
          ("vof_floatingCylinder_3d_p" , #2
           "vof_floatingCylinder_3d_n"),
          ("redist_floatingCylinder_3d_p" ,#3 
           "redist_floatingCylinder_3d_n"),
          ("ls_consrv_floatingCylinder_3d_p" ,#4 
           "ls_consrv_floatingCylinder_3d_n"),
          ("moveMesh_floatingCylinder_3d_p",#5
           "moveMesh_floatingCylinder_3d_n")]

name = "floatingCylinder"

systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False
useOneArchive = False
archiveFlag = ArchiveFlags.EVERY_USER_STEP

if floatingCylinder.nDTout != None:
    tnList = [0.0,floatingCylinder.dt_init] + \
    [floatingCylinder.dt_init+i*(floatingCylinder.T-floatingCylinder.dt_init)/(float(floatingCylinder.nDTout)-1) 
     for i in range(1,floatingCylinder.nDTout)]
else:
    tnList = [0.0,floatingCylinder.dt_init,floatingCylinder.T]
