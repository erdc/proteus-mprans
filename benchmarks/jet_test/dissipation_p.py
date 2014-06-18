from proteus import *
from proteus.default_p import *
from jet_test3D import *
from proteus.mprans import Dissipation

LevelModelType = Dissipation.LevelModel

coefficients = Dissipation.Coefficients(V_model=0,ME_model=4,LS_model=1,RD_model=None,kappa_model=3,
                                        dissipation_model_flag=dissipation_model_flag,
                                        useMetrics=useMetrics,
                                        rho_0=rho_0,nu_0=nu_0,
                                        rho_1=rho_1,nu_1=nu_1,
                                        g=g,
                                        c_mu=0.09,sigma_e=1.0,#default values for c_1,c_2,c_e
                                        sc_uref=dissipation_sc_uref,sc_beta=dissipation_sc_beta)

dissipationInflow = coefficients.c_mu*kInflow**(1.5)/(0.03*inflow_length)
if dissipation_model_flag >= 2:
    dissipationInflow = dissipationInflow/(kInflow+1.0e-12)

def getDBC_dissipation(x,flag):
    if flag == boundaryTags['inflow']:
        return lambda x,t:dissipationInflow

dirichletConditions = {0:getDBC_dissipation}
fluxBoundaryConditions = {0:'outFlow'}

def getAFBC_dissipation(x,flag):
    if flag == boundaryTags['outlet']:
        return None
    if flag in wall_boundaries:
        return lambda x,t: 0.0#no flux
def getDFBC_dissipation(x,flag):
    if flag == boundaryTags['outlet']:
        return lambda x,t: 0.0#outflow
    if flag != boundaryTags['inflow']:
        return lambda x,t: 0.0#outflow or no flux
    

advectiveFluxBoundaryConditions =  {0:getAFBC_dissipation}
diffusiveFluxBoundaryConditions = {0:{0:getDFBC_dissipation}}



class ConstantIC:
    def __init__(self,cval=0.0):
        self.cval=cval
    def uOfXT(self,x,t):
        return self.cval
   
initialConditions  = {0:ConstantIC(cval=dissipationInflow*0.001)}
