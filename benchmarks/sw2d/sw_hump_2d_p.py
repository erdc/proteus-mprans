from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2D
nd=2


L=(2.0,2.0,1.0)
g = 1.0
H0=1.0
HH=2.0
shock=True

#coefficients = ShallowWater(g=g,
#                            nd=nd)
LevelModelType = SW2D.LevelModel
coefficients = SW2D.Coefficients(nu=0.0,g=1.0)
#coefficients = SW2D.Coefficients(nu=0,g=1.0)

class HumpIC:
    def __init__(self,Lx,Ly,r,H,H0,shock=False):
        self.shock=shock
        self.xc = Lx/2.0
        self.yc = Ly/2.0
        self.r=r
        self.H=H
        self.H0=H0
    def uOfXT(self,X,t):
        import math
        x = X[0]
        y = X[1]
        if math.sqrt((x - self.xc)**2 + (y - self.yc)**2) < self.r:
            if self.shock:
                h = self.H
            else:
                h = (self.H*
                     (1.0+cos(math.pi*(x-self.xc)/self.r))*
                     (1.0+cos(math.pi*(y-self.yc)/self.r))
                     +self.H0)
        else:
            h = self.H0
        return h

class ZeroIC:
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:HumpIC(L[0],L[1],0.25*L[0],HH,H0),
                     1:ZeroIC(),
                     2:ZeroIC()}

def getDBC_h(x,flag):
    return None

def getDBC_u(x,flag):
#    return None
   if x[0] in [0.0,L[0]]:
       return lambda x,t: 0.0
   else:
       return None

def getDBC_v(x,flag):
#    return None
   if x[1] in [0.0,L[1]]:
       return lambda x,t: 0.0
   else:
       return None

dirichletConditions = {0:getDBC_h,
                       1:getDBC_u,
                       2:getDBC_v}

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

def getAFBC_h(x,flag):
    return lambda x,t: 0.0

def getAFBC_u(x,flag):
#    return lambda x,t: 0.0
    return None
def getAFBC_v(x,flag):
#    return lambda x,t: 0.0
    return None

advectiveFluxBoundaryConditions =  {0:getAFBC_h,
                                    1:getAFBC_u,
                                    2:getAFBC_v}

def getDFBC_u(x,flag):
    return None
#    return lambda x,t: 0.0

def getDFBC_v(x,flag):
    return None
#    return lambda x,t: 0.0

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u},
                                   2:{2:getDFBC_v}}

T=1.0


