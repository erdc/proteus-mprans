from proteus.default_so import *
import obstacleInTank3d

pnList = [("twp_navier_stokes_obstacleInTank_3d_p","twp_navier_stokes_obstacleInTank_3d_n"),
          ("ls_obstacleInTank_3d_p" ,              "ls_obstacleInTank_3d_n"),
          ("vof_obstacleInTank_3d_p" ,             "vof_obstacleInTank_3d_n"),
          ("redist_obstacleInTank_3d_p" ,          "redist_obstacleInTank_3d_n"),
          ("ls_consrv_obstacleInTank_3d_p" ,       "ls_consrv_obstacleInTank_3d_n")]

    
name = "twp_navier_stokes_obstacleInTank_3d"

if obstacleInTank3d.useBackwardEuler:
    systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinFLCBDFModelStep
if obstacleInTank3d.useFixedStep:
    systemStepControllerType = Sequential_FixedStep_Simple 
needEBQ_GLOBAL = False#True
needEBQ = False#True
useOneArchive = False

tnList = [0.0,obstacleInTank3d.dt_init]+[obstacleInTank3d.dt_init+i*(obstacleInTank3d.T-obstacleInTank3d.dt_init)/float(obstacleInTank3d.nDTout-1) for i in range(1,obstacleInTank3d.nDTout)]
if obstacleInTank3d.useFixedStep:
    tnList = [0.0,obstacleInTank3d.dt_init]+[(i+1)*obstacleInTank3d.dt_fixed for i in range(1,obstacleInTank3d.nDTout)]
