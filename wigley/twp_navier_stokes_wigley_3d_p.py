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
                                   stokes=useStokes,
                                   useRBLES=useRBLES,
		                   useMetrics=useMetrics)
coefficients.waterLevel=waterLevel

def velRamp(t):
    if rampInitialConditions and (t < 0.5*residence_time):
        return (t/(0.5*residence_time))
    else:
        return 1.0

def getDBC_p(x,flag):
    def p(x,t):
        return -coefficients.g[2]*(rho_0*(height - x[2])
                                   -(rho_0-rho_1)*smoothedHeaviside_integral(epsFact_density*he,height-waterLevel)
                                   +(rho_0-rho_1)*smoothedHeaviside_integral(epsFact_density*he,x[2]-waterLevel))
    if flag == boundaryTags['right']:
        return p
    if flag == boundaryTags['front']:
        if openSides: return p
    if flag == boundaryTags['back']:
        if openSides: return p
    if flag == boundaryTags['top']:
        if openTop: return p
    
def getDBC_u(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: Um*velRamp(t)
    if flag == boundaryTags['front']:
        if openSides: return lambda x,t: Um*velRamp(t)
    if flag == boundaryTags['back']:
        if openSides: return lambda x,t: Um*velRamp(t)
    if flag == boundaryTags['top']:
        if openTop: return lambda x,t: Um*velRamp(t)
    if flag == boundaryTags['bottom']:
        if not smoothBottom: return lambda x,t: Um*velRamp(t)
    if flag == boundaryTags['obstacle']:
        if not smoothObstacle:
            if movingDomain:
                return lambda x,t: rc.get_u()
            else:
                return lambda x,t: 0.0

def getDBC_v(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: 0.0
    if flag == boundaryTags['top']:
        if openTop: return lambda x,t: 0.0
    if flag == boundaryTags['bottom']:
        if not smoothBottom: return lambda x,t: 0.0
    if flag == boundaryTags['obstacle']:
        if not smoothObstacle:
            if movingDomain:
                return lambda x,t: rc.get_v()
            else:
                return lambda x,t: 0.0

def getDBC_w(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: 0.0
#     if flag == boundaryTags['front']:
#         if openSides: return lambda x,t: 0.0
#     if flag == boundaryTags['back']:
#         if openSides: return lambda x,t: 0.0
    if flag == boundaryTags['bottom']:
        if not smoothBottom: return lambda x,t: 0.0
    if flag == boundaryTags['obstacle']:
        if not smoothObstacle:
            if movingDomain:
                return lambda x,t: rc.get_w()
            else:
                return lambda x,t: 0.0

dirichletConditions = {0:getDBC_p,
                       1:getDBC_u,
                       2:getDBC_v,
                       3:getDBC_w}

def getAFBC_p(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: -Um*velRamp(t)
#    if flag == boundaryTags['right']:    #IDO
#        return lambda x,t: -Um*velRamp(t)
    if flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0
    if flag == boundaryTags['bottom']:
        return lambda x,t: 0.0
    if flag == boundaryTags['top']:
        if not openTop: return lambda x,t: 0.0
    if flag == boundaryTags['front']:
        if not openSides: return lambda x,t: 0.0
    if flag == boundaryTags['back']:
        if not openSides: return lambda x,t: 0.0
    #processor boundaries in parallel
    if flag == 0:
        return lambda x,t: 0.0

def getAFBC_u(x,flag):
    if flag == boundaryTags['top']:
        if not openTop: return lambda x,t: 0.0
    if flag == boundaryTags['front']:
        if not openSides: return lambda x,t: 0.0
    if flag == boundaryTags['back']:
        if not openSides: return lambda x,t: 0.0
    if flag == boundaryTags['obstacle']:
        if smoothObstacle: return lambda x,t: 0.0
    if flag == boundaryTags['bottom']:
        if smoothBottom: return lambda x,t: 0.0
    #processor boundaries in parallel
    if flag == 0:
        return lambda x,t: 0.0

def getAFBC_v(x,flag):
    if flag == boundaryTags['top']:
        if not openTop: return lambda x,t: 0.0
    if flag == boundaryTags['front']:
        if not openSides: return lambda x,t: 0.0
    if flag == boundaryTags['back']:
        if not openSides: return lambda x,t: 0.0
    if flag == boundaryTags['bottom']:
        if smoothBottom: return lambda x,t: 0.0
    if flag == boundaryTags['obstacle']:
        if smoothObstacle: return lambda x,t: 0.0
    #processor boundaries in parallel
    if flag == 0:
        return lambda x,t: 0.0
    
def getAFBC_w(x,flag):
    if flag == boundaryTags['top']:
        if not openTop: return lambda x,t: 0.0
    if flag == boundaryTags['front']:
        if not openSides: return lambda x,t: 0.0
    if flag == boundaryTags['back']:
        if not openSides: return lambda x,t: 0.0
    if flag == boundaryTags['obstacle']:
        if smoothObstacle: return lambda x,t: 0.0
    if flag == boundaryTags['bottom']:
        if smoothBottom: return lambda x,t: 0.0
    #processor boundaries in parallel
    if flag == 0:
        return lambda x,t: 0.0
    
def getDFBC_u(x,flag):
    if flag == boundaryTags['top']:
        if not openTop: return lambda x,t: 0.0 #no flow
    if flag == boundaryTags['front']:
        if not openSides: return lambda x,t: 0.0 #no flow
    if flag == boundaryTags['back']:
        if not openSides: return lambda x,t: 0.0 #no flow
    if flag == boundaryTags['right']:
        return lambda x,t: 0.0 #outflow
    if flag == boundaryTags['obstacle']:
        if smoothObstacle: return lambda x,t: 0.0
    if flag == boundaryTags['bottom']:
        if smoothBottom: return lambda x,t: 0.0
    #processor boundaries in parallel
    if flag == 0:
        return lambda x,t: 0.0

def getDFBC_v(x,flag):
    if flag == boundaryTags['top']:
        return lambda x,t: 0.0 #no flow or outflow
    if flag == boundaryTags['right']:
        return lambda x,t: 0.0 #outflow
    if flag == boundaryTags['front']:
        return lambda x,t: 0.0 #no flow or outflow
    if flag == boundaryTags['back']:
        return lambda x,t: 0.0 #no flow or outflow
    if flag == boundaryTags['obstacle']:
        if smoothObstacle: return lambda x,t: 0.0
    if flag == boundaryTags['bottom']:
        if smoothBottom: return lambda x,t: 0.0
    #processor boundaries in parallel
    if flag == 0:
        return lambda x,t: 0.0

def getDFBC_w(x,flag):
    #print  boundaryTags
    if flag == boundaryTags['top']:
        return lambda x,t: 0.0 #no flow or outflow
    if flag == boundaryTags['right']:
        return lambda x,t: 0.0 #outflow
    if flag == boundaryTags['front']:
        return lambda x,t: 0.0 #no flow or outflow
        #if not openSides: return lambda x,t: 0.0 #no flow or outflow
    if flag == boundaryTags['back']:
        return lambda x,t: 0.0 #no flow or outflow
        #if not openSides: return lambda x,t: 0.0 #no flow or outflow
    if flag == boundaryTags['obstacle']:
        if smoothObstacle: return lambda x,t: 0.0
    if flag == boundaryTags['bottom']:
        if smoothBottom: return lambda x,t: 0.0
    #processor boundaries in parallel
    if flag == 0:
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
                                   1:{1:getDFBC_u},
                                   2:{2:getDFBC_v},
                                   3:{3:getDFBC_w}}

class Steady_p:
    def uOfXT(self,x,t):
        return -coefficients.g[2]*(rho_0*(height - x[2])
                                   -(rho_0-rho_1)*smoothedHeaviside_integral(epsFact_density*he,height-waterLevel)
                                   +(rho_0-rho_1)*smoothedHeaviside_integral(epsFact_density*he,x[2]-waterLevel))

class Steady_u:
    def uOfXT(self,x,t):
        return velRamp(t)*Um
    
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
