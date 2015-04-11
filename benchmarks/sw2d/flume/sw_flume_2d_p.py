from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2D
from proteus.mprans import SW2DCV
nd=2


L=(4.0,2.0)#needs to be 2d
g = 1.0
nu= 1.0e-8#1.0e-3#0.1
H1=2.0
H0=1.0
Q=2.0e-1
shock=True
domain = Domain.RectangularDomain(L,name="Rectangle")
bt = domain.boundaryTags
#domain.writePoly("domain2d")

LevelModelType = SW2D.LevelModel
coefficients = SW2D.Coefficients(nu=nu,g=g)

class ShockIC:
    def __init__(self,val0=0.0,val1=0.0):
        self.val0 = val0
        self.val1 = val1
    def uOfXT(self,x,t):
        if x[0] <= L[0]*0.5:
            return self.val0
        return self.val1

class ConstIC:
    def __init__(self,val=0.0):
        self.val = val
    def uOfXT(self,x,t):
        return self.val

initialConditions = {0:ShockIC(H0,H0),
                     1:ConstIC(),
                     2:ConstIC()}

def getDBC_h(x,flag):
    #return None
    if flag == bt['right']:
        return lambda x,t: H0
    return None
def getDBC_u(x,flag):
    if flag in [bt['top'],bt['bottom']]:#,bt['left']]:#,bt['right']]:
        return lambda x,t: 0.0
    if flag in [bt['left']]:
        return lambda x,t: Q
    return None
    


def getDBC_v(x,flag):
    if flag in [bt['top'],bt['bottom'],bt['left'],bt['right']]:
        return lambda x,t: 0.0
    return None

dirichletConditions = {0:getDBC_h,
                       1:getDBC_u,
                       2:getDBC_v}

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

def getAFBC_h(x,flag):
    #if flag in [bt['top'],bt['bottom'],bt['left'],bt['right']]:
    #    return lambda x,t: 0.0
    if flag not in [bt['left'],bt['right']]:
        return lambda x,t: 0.0
    if flag == bt['left']:
        return lambda x,t: -Q

#not right yet
def getAFBC_u(x,flag):
    return None
    return lambda x,t: 0.0

def getAFBC_v(x,flag):
    return None
    return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_h,
                                    1:getAFBC_u,
                                    2:getAFBC_v}

#ignored right now
def getDFBC_u(x,flag):
    return None

def getDFBC_v(x,flag):
    return None

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u},
                                   2:{2:getDFBC_v}}

T=5.0


