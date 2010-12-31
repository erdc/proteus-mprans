run with:
parun wigley_so.py
files:
README.txt : this file
ls_consrv_wigley_3d_n.py : level set conservation numerics
ls_consrv_wigley_3d_p.py : level set conservation physics
ls_wigley_3d_n.py : level set advection numerics
ls_wigley_3d_p.py : level set advection physics
redist_wigley_3d_n.py : eikonal (redistancing) numerics
redist_wigley_3d_p.py : eikonal (redistancing) phyics
twp_navier_stokes_wigley_3d_n.py : flow numerics
twp_navier_stokes_wigley_3d_p.py : flow physics
vof_wigley_3d_n.py :  volume fraction numerics
vof_wigley_3d_p.py :  volume fraction physics
wigley.py : "global" input module shared by all p and n files
wigley.ranger.sge : batch script for ranger at TACC
wigley3d.diamond.pbs : batch script for diamond at ERDC DSRC
wigley3d.pbs : batch script for jade at ERDC DSRC
wigley3dDomain.py : domain geometry
wigley_so.py : split operator specification for entire system
