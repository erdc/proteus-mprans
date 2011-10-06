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
from proteus.mprans import RBLES

LevelModelType = RBLES.LevelModel
coefficients = RBLES.Coefficients(epsFact=epsFact_viscosity,
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
coefficients.waterLevel=waterLevel

def velRamp(t):
    return 1.0
    #if rampInitialConditions and (t < 0.5*residence_time):
    #    return (t/(0.5*residence_time))
    #else:
    #    return 1.0
   
#----------------------------------------------------------
#  Weak dirichlet
#----------------------------------------------------------    
def getDBC_p(x,flag):
    pass
  
def getDBC_u(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: Um*velRamp(t)
    if flag == boundaryTags['bottom']:
        return lambda x,t: Um*velRamp(t)
    if flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0
	
def getDBC_v(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: 0.0
    if flag == boundaryTags['bottom']:
        return lambda x,t: 0.0
    if flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0
	
def getDBC_w(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: 0.0
    if flag == boundaryTags['bottom']:
        return lambda x,t: 0.0
    if flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0
	
dirichletConditions = {0:getDBC_p,
                       1:getDBC_u,
                       2:getDBC_v,
                       3:getDBC_w}

#----------------------------------------------------------
#  Do nothing, convective fallback values...
#----------------------------------------------------------   
def getAFBC_p(x,flag):
    pass 

def getAFBC_u(x,flag):
    if flag == boundaryTags['top']:
        return lambda x,t: Um*velRamp(t)
    if flag == boundaryTags['front']:
        return lambda x,t: Um*velRamp(t)
    if flag == boundaryTags['back']:
        return lambda x,t: Um*velRamp(t)
    if flag == boundaryTags['right']:
        return lambda x,t: Um*velRamp(t)

def getAFBC_v(x,flag):
    if flag == boundaryTags['top']:
        return lambda x,t: 0.0
    if flag == boundaryTags['front']:
        return lambda x,t: 0.0
    if flag == boundaryTags['back']:
        return lambda x,t: 0.0
    if flag == boundaryTags['right']:
        return lambda x,t: 0.0

def getAFBC_w(x,flag):
    if flag == boundaryTags['top']:
        return lambda x,t: 0.0
    if flag == boundaryTags['front']:
        return lambda x,t: 0.0
    if flag == boundaryTags['back']:
        return lambda x,t: 0.0
    if flag == boundaryTags['right']:
        return lambda x,t: 0.0	

def getDFBC(x,flag):
        return lambda x,t: 0.0

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
#  Initial condtionss
#----------------------------------------------------------
class Steady_p:
    def uOfXT(self,x,t):
        return 0.0  

class Steady_u:
    def uOfXT(self,x,t):
        return Um*velRamp(t)
    
class Steady_v:
    def uOfXT(self,x,t):
        return 0.0

class Steady_w:
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:Steady_p(),
                     1:Steady_u(),
                     2:Steady_v(),
                     3:Steady_w()}
