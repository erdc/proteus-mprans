from proteus import *
from proteus.default_p import *
from dtmb import *
from proteus.mprans import Epsilon

LevelModelType = Epsilon.LevelModel

if useOnlyVF:
    coefficients = Epsilon.Coefficients(V_model=0,ME_model=3,LS_model=None,RD_model=None,kappa_model=2,
                                        useMetrics=useMetrics,
                                        rho_0=rho_0,nu_0=nu_0,
                                        rho_1=rho_1,nu_1=nu_1,
                                        g=g,
                                        c_mu=0.09,sigma_e=1.0,#default values for c_1,c_2,c_e
                                        sc_uref=epsilon_sc_uref,sc_beta=epsilon_sc_beta)
    
else:
    coefficients = Epsilon.Coefficients(V_model=0,ME_model=6,LS_model=2,RD_model=3,kappa_model=5,
                                        useMetrics=useMetrics,
                                        rho_0=rho_0,nu_0=nu_0,
                                        rho_1=rho_1,nu_1=nu_1,
                                        g=g,
                                        c_mu=0.09,sigma_e=1.0,#default values for c_1,c_2,c_e
                                        sc_uref=epsilon_sc_uref,sc_beta=epsilon_sc_beta)

epsilonInflow = coefficients.c_mu*kInflow**(1.5)/(0.03*L[2])

def getDBC_epsilon(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t:epsilonInflow

dirichletConditions = {0:getDBC_epsilon}
fluxBoundaryConditions = {0:'outFlow'}

def getAFBC_epsilon(x,flag):
    if flag == boundaryTags['right']:
        return None
    if flag != boundaryTags['left']:
        return lambda x,t: 0.0
def getDFBC_epsilon(x,flag):
    pass
#    if flag == boundaryTags['right']:
#        return lambda x,t: 0.0
#    if flag != boundaryTags['left']:
#        return lambda x,t: 0.0
    

advectiveFluxBoundaryConditions =  {0:getAFBC_epsilon}
diffusiveFluxBoundaryConditions = {0:{0:getDFBC_epsilon}}



class ConstantIC:
    def __init__(self,cval=0.0):
        self.cval=cval
    def uOfXT(self,x,t):
        return self.cval
   
initialConditions  = {0:ConstantIC(cval=epsilonInflow*0.001)}
