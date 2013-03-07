from proteus import *
from proteus.default_p import *
from dtmb import *
from proteus.mprans import Kappa

LevelModelType = Kappa.LevelModel
if useOnlyVF:
    coefficients = Kappa.Coefficients(V_model=0,ME_model=2,LS_model=None,RD_model=None,dissipation_model=3,
                                      useMetrics=useMetrics,
                                      rho_0=rho_0,nu_0=nu_0,
                                      rho_1=rho_1,nu_1=nu_1,
                                      g=g,
                                      c_mu=0.09,sigma_k=1.0,
                                      sc_uref=kappa_sc_uref,sc_beta=kappa_sc_beta)
else:
    coefficients = Kappa.Coefficients(V_model=0,ME_model=5,LS_model=2,RD_model=3,dissipation_model=6,
                                      useMetrics=useMetrics,
                                      rho_0=rho_0,nu_0=nu_0,
                                      rho_1=rho_1,nu_1=nu_1,
                                      g=g,
                                      c_mu=0.09,sigma_k=1.0,
                                      sc_uref=kappa_sc_uref,sc_beta=kappa_sc_beta)


def getDBC_k(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t:kInflow

dirichletConditions = {0:getDBC_k}
fluxBoundaryConditions = {0:'outFlow'}

def getAFBC_k(x,flag):
    if flag == boundaryTags['right']:
        return None
    if flag != boundaryTags['left']:
        return lambda x,t: 0.0
def getDFBC_k(x,flag):
    pass
#    if flag == boundaryTags['right']:
#        return lambda x,t: 0.0
#    if flag != boundaryTags['left']:
#        return lambda x,t: 0.0
    

advectiveFluxBoundaryConditions =  {0:getAFBC_k}
diffusiveFluxBoundaryConditions = {0:{0:getDFBC_k}}



class ConstantIC:
    def __init__(self,cval=0.0):
        self.cval=cval
    def uOfXT(self,x,t):
        return self.cval
   
initialConditions  = {0:ConstantIC(cval=kInflow*0.001)}
