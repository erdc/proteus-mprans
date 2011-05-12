run with:
parun obstacleInTank3d_so.py 

files:
boxInTank3dDomain : domain geometry
README.txt : this file
ls_consrv_obstacleInTank_3d_n.py : level set conservation numerics
ls_consrv_obstacleInTank_3d_p.py : level set conservation physics
ls_obstacleInTank_3d_n.py : level set numerics
ls_obstacleInTank_3d_p.py : level set physics
obstacle.lonestar.lsf : batch script for lonestar at TACC
obstacleInTank3d.pbs : batch script for jade at ERDC DSRC
obstacleInTank3d.py : "global" inputs
obstacleInTank3d_so.py : split operator specification for system
redist_obstacleInTank_3d_n.py : eikonal (redistancing equation) numerics
redist_obstacleInTank_3d_p.py : eikonal (redistancing equation) physics
twp_navier_stokes_obstacleInTank_3d_n.py : flow numerics
twp_navier_stokes_obstacleInTank_3d_p.py : flow physics
vof_obstacleInTank_3d_n.py : volume fraction numerics
vof_obstacleInTank_3d_p.py : volume fraction physics
