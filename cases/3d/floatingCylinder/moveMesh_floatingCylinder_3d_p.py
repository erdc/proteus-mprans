from proteus import *
from proteus.default_p import *
from floatingCylinder import *
from proteus.mprans import MoveMesh

initialConditions = None

analyticalSolution = {}

LevelModelType = MoveMesh.LevelModel

nMediaTypes=1
smTypes      = numpy.zeros((nMediaTypes,2),'d')
smFlags      = numpy.zeros((nMediaTypes,),'i')

E=he**3/6.0#1.0e5 #kN/m^2 Young's modulus
nu=0.3  #- Poisson's ratio

smTypes[0,0] = E
smTypes[0,1] = nu
coefficients = MoveMesh.Coefficients(hullMass=0.0,
                                     hullCG=[0.0,0.0,0.0],
                                     hullInertia=[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]],
                                     modelType_block=smFlags,modelParams_block=smTypes,meIndex=5)
#coefficients = MoveMesh.Coefficients(E=10.0,nu=0.3,g=g,nd=nd,moveMesh=movingDomain)

fixedBoundary =  [boundaryTags['bottom'],boundaryTags['top'],boundaryTags['front'],boundaryTags['back'],boundaryTags['upstream'],boundaryTags['downstream']]

class TranslatingObstacle(AuxiliaryVariables.AV_base):
    def __init__(self):
        self.current_center=cylinder_center
        self.r = 0.0*cylinder_radius
        self.omega=1.0#1.0/1.0
        self.cx = cylinder_center[0]
        self.cy = cylinder_center[1]
        self.cz = cylinder_center[2]
        self.h=[0.0,0.0,0.0]
    def attachModel(self,model,ar):
        self.model=model
        return self
    def attachAuxiliaryVariables(self,avDict):
        self.object = avDict['twp_navier_stokes_floatingCylinder_3d_p'][0]
    def center(self,t):
        return (self.cx+self.r*cos(self.omega*2.0*math.pi*t+math.pi),
                self.cy,
                self.cz+self.r*sin(self.omega*2.0*math.pi*t +math.pi))
    def hx(self,t):
        h = tro.center(t)[0]-tro.current_center[0]
        return h
    def hy(self,t):
        h = tro.center(t)[1]-tro.current_center[1]
        return h
    def hz(self,t):
        h = tro.center(t)[2]-tro.current_center[2]
        return h
    def calculate(self):
        self.h[0] = self.hx(self.model.stepController.t_model)
        self.h[1] = self.hy(self.model.stepController.t_model)
        self.h[2] = self.hz(self.model.stepController.t_model)
        print "displacement------------------------------------",self.h[0],self.h[1],self.h[2]
        self.object.last_velocity = (self.h[0]/self.model.stepController.dt_model,
                                     self.h[1]/self.model.stepController.dt_model,
                                     self.h[2]/self.model.stepController.dt_model)
        self.current_center=self.center(self.model.stepController.t_model_last)

class TranslatingFreeObstacle(AuxiliaryVariables.AV_base):
    def __init__(self):
        self.current_center=cylinder_center
        self.r = 0.5*cylinder_radius
        self.omega=1.0/1.0
        self.cx = cylinder_center[0]
        self.cy = cylinder_center[1]
        self.cz = cylinder_center[2]
        self.body=None
        self.h=(0.0,0.0,0.0)
        self.object = None
    def attachModel(self,model,ar):
        self.model=model
        return self
    def attachAuxiliaryVariables(self,avDict):
        self.object = avDict['twp_navier_stokes_floatingCylinder_3d_p'][0]
    def hx(self,t):
        if self.object == None:
            return 0.0
        else:
            return self.object.h[0]
    def hy(self,t):
        if self.object == None:
            return 0.0
        else:
            return self.object.h[1]
    def hz(self,t):
        if self.object == None:
            return 0.0
        else:
            return self.object.h[2]
    def calculate(self):
        self.h=self.object.h#(self.object.position[0]-self.object.last_position[0],
        print "displacement------------------------------------",self.h[0],self.h[1],self.h[2]

tro = TranslatingFreeObstacle()
tro = TranslatingObstacle()

def getDBC_hx(x,flag):
    if flag in [boundaryTags['upstream'],boundaryTags['downstream']]:
        return lambda x,t: 0.0
    if flag == boundaryTags['obstacle']:
        return lambda x,t: tro.hx(t)
    
def getDBC_hy(x,flag):
    if flag in [boundaryTags['front'],boundaryTags['back']]:
        return lambda x,t: 0.0
    if flag == boundaryTags['obstacle']:
        return lambda x,t: tro.hy(t)
    
def getDBC_hz(x,flag):
    if flag in [boundaryTags['top'],boundaryTags['bottom']]:
        return lambda x,t: 0.0
    if flag == boundaryTags['obstacle']:
        return lambda x,t: tro.hz(t)

dirichletConditions = {0:getDBC_hx,
                       1:getDBC_hy,
                       2:getDBC_hz}

fluxBoundaryConditions = {0:'noFlow',
                          1:'noFlow',
                          2:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{},
                                   2:{}}

def stress_u(x,flag):
    if flag not in [boundaryTags['upstream'],
                    boundaryTags['downstream']]:
        return 0.0
    
def stress_v(x,flag):
    if flag not in [boundaryTags['front'],
                    boundaryTags['back']]:
        return 0.0

def stress_w(x,flag):
    if flag not in [boundaryTags['top'],
                    boundaryTags['bottom']]:
        return 0.0

stressFluxBoundaryConditions = {0:stress_u,
                                1:stress_v,
                                2:stress_w}
