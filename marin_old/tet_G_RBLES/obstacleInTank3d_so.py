from proteus.default_so import *
import obstacleInTank3d

pnList = [("twp_navier_stokes_obstacleInTank_3d_p","twp_navier_stokes_obstacleInTank_3d_n"),
          ("ls_obstacleInTank_3d_p" ,              "ls_obstacleInTank_3d_n"),
          ("vof_obstacleInTank_3d_p" ,             "vof_obstacleInTank_3d_n"),
          ("redist_obstacleInTank_3d_p" ,          "redist_obstacleInTank_3d_n"),
          ("ls_consrv_obstacleInTank_3d_p" ,       "ls_consrv_obstacleInTank_3d_n")]

    
name = "twp_navier_stokes_obstacleInTank_3d"

systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False
useOneArchive = False

tnList = [0.0,obstacleInTank3d.dt_init]+[i*0.05 for i in range(1,120)]

