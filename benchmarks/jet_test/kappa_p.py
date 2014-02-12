from proteus import *
from proteus.default_p import *
from jet_test2D import *
from proteus.mprans import Kappa

LevelModelType = Kappa.LevelModel

coefficients = Kappa.Coefficients(V_model=0,ME_model=3,LS_model=1,RD_model=None,dissipation_model=4,
                                  dissipation_model_flag=dissipation_model_flag,
                                  useMetrics=useMetrics,
                                  rho_0=rho_0,nu_0=nu_0,
                                  rho_1=rho_1,nu_1=nu_1,
                                  g=g,
                                  c_mu=0.09,sigma_k=1.0,
                                  sc_uref=kappa_sc_uref,sc_beta=kappa_sc_beta)


def getDBC_k(x,flag):
    if flag == boundaryTags['inflow']:
        return lambda x,t:kInflow
    if flag == boundaryTags['wall']:
        return lambda x,t:0.0
dirichletConditions = {0:getDBC_k}
fluxBoundaryConditions = {0:'outFlow'}

def getAFBC_k(x,flag):
    if flag == boundaryTags['outflow']:
        return None
    if flag != boundaryTags['inflow']:
        return None#outflow

def getDFBC_k(x,flag):
    if flag == boundaryTags['outflow']:
        return lambda x,t: 0.0#outflow
    if flag != boundaryTags['inflow']:
        return lambda x,t: 0.0#outflow or no flux
    

advectiveFluxBoundaryConditions =  {0:getAFBC_k}
diffusiveFluxBoundaryConditions = {0:{0:getDFBC_k}}



class ConstantIC:
    def __init__(self,cval=0.0):
        self.cval=cval
    def uOfXT(self,x,t):
        return self.cval
   
initialConditions  = {0:ConstantIC(cval=kInflow*0.001)}
