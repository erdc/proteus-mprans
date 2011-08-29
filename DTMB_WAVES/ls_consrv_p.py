from proteus import *
from proteus.default_p import *
from DTMB import *
from proteus.mprans import MCorr

LevelModelType = MCorr.LevelModel
coefficients = MCorr.Coefficients(applyCorrection=applyCorrection,
                                    epsFactHeaviside=epsFact_consrv_heaviside,
                                    epsFactDirac=epsFact_consrv_dirac,
                                    epsFactDiffusion=epsFact_consrv_diffusion,
                                    LSModel_index=1,
                                    V_model=0,
                                    me_model=4,
                                    VOFModel_index=2,
                                    nd=nd,
                                    checkMass=checkMass,
		                    useMetrics=useMetrics)

class zero_phi:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return vof_wave(x,t)

analyticalSolutions = None

def getDBC_cnsrv(x,flag):
    pass

dirichletConditions = {0:getDBC_cnsrv}
#bubble rise
initialConditions  = {0:zero_phi()}

fluxBoundaryConditions = {0:'outFlow'}


def getAFBC_cnsrv(x,flag):
    pass

def getDFBC_cnsrv(x,flag):
    return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_cnsrv}

diffusiveFluxBoundaryConditions = {0:{0:getDFBC_cnsrv}}
