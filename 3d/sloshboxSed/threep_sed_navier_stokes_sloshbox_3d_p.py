from proteus import *
from proteus.default_p import *
from sloshbox3d import *
from proteus.mprans import RANS3PSed

useOpt=True#False
if useOpt:
    LevelModelType = RANS3PSed.LevelModel
coefficients = RANS3PSed.Coefficients(epsFact=epsFact_viscosity,
                                      sigma=0.0,
                                      rho_0 = rho_0,
                                      nu_0 = nu_0,
                                      rho_1 = rho_1,
                                      nu_1 = nu_1,
                                      g=g,
                                      nd=nd,
                                      RANS_model_ID=0,
                                      LS_model=1,
                                      epsFact_density=epsFact_density,
                                      stokes=useStokes)
C0=0.0001

def getDBC_C_sloshbox(x,flag):
   if regularGrid:
        if x[2] > L[2] - 1.0e-8:
            return lambda x,t: C0
   else:
       if flag == bt['top']:
           return lambda x,t: C0
           
def getDBC_us_sloshbox(x,flag):
    if regularGrid:
        if (x[0] == 0.0 or
            x[0] == L[0] or
            x[2] == 0.0):
            return lambda x,t: 0.0
    else:
        if (flag == bt['left'] or
            flag == bt['right'] or
            flag == bt['bottom']):
            return lambda x,t: 0.0

def getDBC_vs_sloshbox(x,flag):
    if regularGrid:
        if (x[0] == 0.0 or
            x[0] == L[0] or
            x[2] == 0.0):
            return lambda x,t: 0.0
    else:
        if (flag == bt['left'] or
            flag == bt['right'] or
            flag == bt['bottom']):
            return lambda x,t: 0.0

def getDBC_ws_sloshbox(x,flag):
    if regularGrid:
        if (x[0] == 0.0 or
            x[0] == L[0] or
            x[2] == 0.0):
            return lambda x,t: 0.0
    else:
        if (flag == bt['left'] or
            flag == bt['right'] or
            flag == bt['bottom']):
            return lambda x,t: 0.0

dirichletConditions = {0:getDBC_C_sloshbox,
                       1:getDBC_us_sloshbox,
                       2:getDBC_vs_sloshbox,
                       3:getDBC_ws_sloshbox}

def getAFBC_C_sloshbox(x,flag):
    if regularGrid:
        if x[2] > L[2] - 1.0e-8:
            return None
        else:
            return lambda x,t: 0.0
    else:
        if(flag == bt['top']):
            return None
        else:
            return lambda x,t: 0.0

def getAFBC_us_sloshbox(x,flag):
    if regularGrid:
        if (x[0] == 0.0 or
            x[0] == L[0] or
            x[2] == 0.0):
            return None
        else:
            return lambda x,t: 0.0
    else:
        if (flag == bt['left'] or
            flag == bt['right'] or
            flag == bt['bottom']):
            return None
        else:
            return lambda x,t: 0.0

def getAFBC_vs_sloshbox(x,flag):
    if regularGrid:
        if (x[0] == 0.0 or
            x[0] == L[0] or
            x[2] == 0.0):
            return None
        else:
            return lambda x,t: 0.0
    else:
        if (flag == bt['left'] or
            flag == bt['right'] or
            flag == bt['bottom']):
            return None
        else:
            return lambda x,t: 0.0

def getAFBC_ws_sloshbox(x,flag):
    if regularGrid:
        if (x[0] == 0.0 or
            x[0] == L[0] or
            x[2] == 0.0):
            return None
        else:
            return lambda x,t: 0.0
    else:
        if (flag == bt['left'] or
            flag == bt['right'] or
            flag == bt['bottom']):
            return None
        else:
            return lambda x,t: 0.0

def getDFBC_us_sloshbox(x,flag):
    if regularGrid:
        if (x[0] == 0.0 or
            x[0] == L[0] or
            x[2] == 0.0):
            return None
        else:
            return lambda x,t: 0.0
    else:
        if (flag == bt['left'] or
            flag == bt['right'] or
            flag == bt['bottom']):
            return None
        else:
            return lambda x,t: 0.0
    
def getDFBC_vs_sloshbox(x,flag):
    if regularGrid:
        if (x[0] == 0.0 or
            x[0] == L[0] or
            x[2] == 0.0):
            return None
        else:
            return lambda x,t: 0.0
    else:
        if (flag == bt['left'] or
            flag == bt['right'] or
            flag == bt['bottom']):
            return None
        else:
            return lambda x,t: 0.0

def getDFBC_ws_sloshbox(x,flag):
    if regularGrid:
        if (x[0] == 0.0 or
            x[0] == L[0] or
            x[2] == 0.0):
            return None
        else:
            return lambda x,t: 0.0
    else:
        if (flag == bt['left'] or
            flag == bt['right'] or
            flag == bt['bottom']):
            return None
        else:
            return lambda x,t: 0.0

fluxBoundaryConditions = {0:'mixedFlow',
                          1:'mixedFlow',
                          2:'mixedFlow',
                          3:'mixedFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_C_sloshbox,
                                    1:getAFBC_us_sloshbox,
                                    2:getAFBC_vs_sloshbox,
                                    3:getAFBC_ws_sloshbox}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_us_sloshbox},
                                   2:{2:getDFBC_vs_sloshbox},
                                   3:{3:getDFBC_ws_sloshbox}}

class Constant_C:
    def __init__(self,C0):
        self.C0=C0
    def uOfXT(self,x,t):
        return self.C0

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:Constant_C(0.0001),
                     1:AtRest(),
                     2:AtRest(),
                     3:AtRest()}
