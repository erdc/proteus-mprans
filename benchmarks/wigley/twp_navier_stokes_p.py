"""
Incompressible Navier-Stokes flow around a square obstacle in 2D.
"""
from proteus import *
from proteus.default_p import *
import sys
from math import *
from wigley import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
from proteus.mprans import RANS2P

LevelModelType = RANS2P.LevelModel
coefficients = RANS2P.Coefficients(epsFact=epsFact_viscosity,
                                   sigma=sigma_01,
                                   rho_0=rho_0,
                                   nu_0=nu_0,
                                   rho_1=rho_1,
                                   nu_1=nu_1,
                                   g=g,
                                   nd=nd,
                                   LS_model=1,
                                   epsFact_density=epsFact_density,
                                   stokes=False,
                                   useRBLES=useRBLES,
		                   useMetrics=useMetrics)
  
#----------------------------------------------------------
#  Weak dirichlet
#----------------------------------------------------------    
def getDBC_p(x,flag):
    if flag == boundaryTags['right']:
      return hs_pres
   # if flag == boundaryTags['top']:
   #   return hs_pres
    
def getDBC_u(x,flag):
    if flag == boundaryTags['left']:
        return u_wave
    if flag == boundaryTags['obstacle']:
        return noslip

def getDBC_v(x,flag):
    if flag == boundaryTags['left']:
        return v_wave
    if flag == boundaryTags['obstacle']:
        return noslip
	
def getDBC_w(x,flag):
    if flag == boundaryTags['left']:
        return w_wave
    if flag == boundaryTags['obstacle']:
        return noslip

dirichletConditions = {0:getDBC_p,
                       1:getDBC_u,
                       2:getDBC_v,
                       3:getDBC_w}

#----------------------------------------------------------
#  Do nothing, convective fallback values...
#----------------------------------------------------------   
def getAFBC_p(x,flag):
    if [ flag == boundaryTags['front']  or
         flag == boundaryTags['back']   or
         flag == boundaryTags['top']    or
         flag == boundaryTags['bottom'] or
         flag == boundaryTags['obstacle'] ]:
        return noflow
    elif flag == boundaryTags['left']:
        return -u_wave
    elif flag == 0:
        return noflow
    	
	
def getAFBC_u(x,flag):
    if [ flag == boundaryTags['front']  or
         flag == boundaryTags['back']   or
         flag == boundaryTags['top']    or
         flag == boundaryTags['bottom'] ]:
        return noflow	
    elif flag == 0:
        return noflow
			
def getAFBC_v(x,flag):
    if [ flag == boundaryTags['front']  or
         flag == boundaryTags['back']   or
         flag == boundaryTags['top']    or
         flag == boundaryTags['bottom'] ]:
        return noflow	
    elif flag == 0:
        return noflow
			
def getAFBC_w(x,flag):
    if [ flag == boundaryTags['front']  or
         flag == boundaryTags['back']   or
         flag == boundaryTags['top']    or
         flag == boundaryTags['bottom'] ]:
        return noflow	
    elif flag == 0:
        return noflow
		
def getDFBC(x,flag):
    if [ flag == boundaryTags['front']  or
         flag == boundaryTags['back']   or
         flag == boundaryTags['top']    or
         flag == boundaryTags['bottom'] or
         flag == boundaryTags['right'] ]:
        return lambda x,t: 0.0
    elif flag == 0:
        return noflow
	
fluxBoundaryConditions = {0:'mixedFlow',
                          1:'mixedFlow',
                          2:'mixedFlow',
                          3:'mixedFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p, 
                                    1:getAFBC_u,
                                    2:getAFBC_v,
                                    3:getAFBC_w}
				    
diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC},
                                   2:{2:getDFBC},
                                   3:{3:getDFBC}}

#----------------------------------------------------------
#  Initial condtions
#----------------------------------------------------------
class Steady_p:
    def uOfXT(self,x,t):
        return hs_pres(x,t)  

class Steady_u:
    def uOfXT(self,x,t):
        return u_wave(x,t)
    
class Steady_v:
    def uOfXT(self,x,t):
        return v_wave(x,t)

class Steady_w:
    def uOfXT(self,x,t):
        return w_wave(x,t)

initialConditions = {0:Steady_p(),
                     1:Steady_u(),
                     2:Steady_v(),
                     3:Steady_w()}
