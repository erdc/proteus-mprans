run with:
parun twp_navier_stokes_sloshbox_3d_so.py 
files:
README.txt : this file
ls_consrv_sloshbox_3d_c0p1_n.py : level set conservation numerics
ls_consrv_sloshbox_3d_p.py : level set conservation physics
ls_sloshbox_3d_c0p1_n.py : level set advection numerics
ls_sloshbox_3d_p.py : level set advection physics
redist_sloshbox_3d_c0p1_n.py : eikonal (redistancing) numerics
redist_sloshbox_3d_p.py : eikonal (redistancing) physics
sloshbox3d.lonestar.lsf : batch script for lonestar at TACC
sloshbox3d.pbs : batch script for jade at ERDC DSRC
sloshbox3d.py : "global" inputs shared  among all p and n modules
sloshbox3d.ranger.sge : batch script for ragner at TACC
twp_navier_stokes_sloshbox_3d_c0p1c0p1_n.py : flow numerics
twp_navier_stokes_sloshbox_3d_p.py : flow physics
twp_navier_stokes_sloshbox_3d_so.py : split operator specification for system
vof_sloshbox_3d_c0p1_n.py : volume fraction numerics
vof_sloshbox_3d_p.py : volume fraction physics
