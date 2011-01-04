from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from vortex import *
from proteus.mprans import VOF
import ls_vortex_3d_p
name=soname+"_vof"

"""
The non-conservative level set description of a bubble in a two-phase flow
"""

LevelModelType = VOF.LevelModel

##\ingroup test
#\file vof_vortex_2d_p.py
#
# \todo finish vof_vortex_2d_p.py

if applyRedistancing:
    coefficients = VOF.Coefficients(LS_model=0,V_model=0,RD_model=1,ME_model=2,checkMass=checkMass,
                                   epsFact=epsFact_vof)
else:
    coefficients = VOF.Coefficients(LS_model=0,V_model=0,RD_model=-1,ME_model=1,checkMass=checkMass,
                                   epsFact=epsFact_vof)

def Heaviside(phi):
    if phi > 0:
        return 1.0
    elif phi < 0:
        return 0.0
    else:
        return 0.5

class H_of_phi:
    def __init__(self,phiSol):
        self.phiSol = phiSol
    def uOfXT(self,X,T):
        return smoothedHeaviside(epsFactHeaviside*he,self.phiSol.uOfXT(X,T))#Heaviside(dBubble)
    def uOfX(self,X):
        return self.uOfXT(X,0)
    #end
#end Vortex_phi

analyticalSolution = {0:H_of_phi(ls_vortex_3d_p.analyticalSolution[0])}

def getDBC(x,flag):
    pass

dirichletConditions = {0:getDBC}

initialConditions  = analyticalSolution

fluxBoundaryConditions = {0:'outFlow'}

#cek made no flux since v.n = 0 for this v
def getAFBC(x,flag):
   return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC}

diffusiveFluxBoundaryConditions = {0:{}}
