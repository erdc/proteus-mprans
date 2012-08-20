from proteus import *
from proteus.default_p import *
from wavetank import *
from proteus.mprans import NCLS

LevelModelType = NCLS.LevelModel

coefficients = NCLS.Coefficients(V_model=0,RD_model=3,ME_model=1,
                                 checkMass=False, useMetrics=useMetrics,
                                 epsFact=epsFact_consrv_heaviside,sc_uref=ls_sc_uref,sc_beta=ls_sc_beta)
 
def getDBC_ls(x,flag):
    if x[0] < 1.0e-8:
        return wavePhi
    if x[0] > L[0] - 1.0e-8:
        return lambda x,t: x[2] - outflowHeight


dirichletConditions = {0:getDBC_ls}

advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{}}

class PHI_IC:       
    def uOfXT(self,x,t):
        return wavePhi(x,t)
    
initialConditions  = {0:PHI_IC()}
