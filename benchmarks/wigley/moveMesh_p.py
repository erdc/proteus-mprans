from proteus import *
from proteus.default_p import *
from wigley  import *
from proteus.mprans import MoveMesh

initialConditions = None

analyticalSolution = {}

nMediaTypes=1
smTypes      = numpy.zeros((nMediaTypes,2),'d')
smFlags      = numpy.zeros((nMediaTypes,),'i')

smTypes[0,0] = 1.0    ##E
smTypes[0,1] = 0.3    ##nu

LevelModelType = MoveMesh.LevelModel
coefficients = MoveMesh.Coefficients(hullMass=hull_mass,    
				     hullCG=hull_cg,      
				     hullInertia=hull_inertia,             	
				     linConstraints=RBR_linCons,  
				     angConstraints=RBR_angCons,  
				     V_model=0,modelType_block=smFlags,
				     modelParams_block=smTypes,meIndex=5)


def getDBC_hx(x,flag):
    if flag in [boundaryTags['top'],boundaryTags['bottom'],boundaryTags['left'],boundaryTags['right'],boundaryTags['front'],boundaryTags['back']]:
        return lambda x,t: 0.0
    if flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0
    
def getDBC_hy(x,flag):
    if flag in [boundaryTags['top'],boundaryTags['bottom'],boundaryTags['left'],boundaryTags['right'],boundaryTags['front'],boundaryTags['back']]:
        return lambda x,t: 0.0
    if flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0
    
def getDBC_hz(x,flag):
    if flag in [boundaryTags['top'],boundaryTags['bottom'],boundaryTags['left'],boundaryTags['right'],boundaryTags['front'],boundaryTags['back']]:
        return lambda x,t: 0.0
    elif flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0
	print "obstacle ",flag


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
    if flag not in [boundaryTags['left'],
                    boundaryTags['right']]:
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
