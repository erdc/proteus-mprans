from proteus import *
from proteus.default_p import *
from math import *
from vortex import *
from proteus.mprans import NCLS
#import Profiling

LevelModelType = NCLS.LevelModel
logEvent = Profiling.logEvent
name=soname+"_ls"

nd=3

## \page Tests Test Problems 
# \ref ls_vortex_2d_p.py "Linear advection of a circular level set function in an oscillating vortex velocity field"
# 

##\ingroup test
# \file la_vortex_2d_p.py
# @{
#  \brief Conservative linear advection of a circle signed distance function
#  in a oscillating vortex velocity field.
#  
# \f{eqnarray*}
# \phi_t + \nabla \cdot (\vec u \phi) &=& 0 \\ 
# \Omega &=& [0,1] \times [0,1] \\
#  u^{x} &=& \cos(\pi t/8)\sin(2\pi y)\sin^2(\pi x) \\
#  u^{y} &=& -\cos(\pi t/8)\sin(2\pi x)\sin^{2}(\pi y) \\
#  \phi^{0} &=& \left(x-\frac{1}{2}\right)^2 + \left(y-\frac{3}{4}\right)^2 - 0.15^2
# \f}
# The solution should return to the initial condition at \f$T=8\f$.
# Outflow boundaries are applied on \f$\partial \Omega\f$.
# 
#
# \image html  save_la_vortex_2d_dgp2_exact.jpg "exact solution, T=8.0"
# \image latex save_la_vortex_2d_dgp2_exact.eps "exact solution, T=8.0"
# \image html  save_la_vortex_2d_dgp2_phi.jpg "RKDG P^2 solution, Cr=0.1, L^2 error= 7.84e-3"
# \image latex save_la_vortex_2d_dgp2_phi.eps "RKDG $P^2$ solution, Cr=0.1, $L^2$ error= 7.84e-3"
#

class OscillatingVortex3D:
    #cek changed to put sphere inside arbitrary box with dimensions in L
    def __init__(self,L):
        self.radius = 0.15*L[0]
        self.xc=0.5*L[0]
        self.yc=0.5*L[1]
        self.zc=0.5*L[2]
    def uOfXT(self,x,t):
        return self.radius - math.sqrt((x[0]-self.xc)**2 + (x[1]-self.yc)**2 + (x[2]-self.zc)**2)
class OscillatingVortex3Dcylinder:
    #cek changed to put sphere inside arbitrary box with dimensions in L
    def __init__(self,L=[1.0,
                         1.0],
                 center=[0.5,
                         0.5],
                 radius=0.45,
                 slotwidth=0.45/3.0,#3/20
                 slotlength=0.45):#45/100 = 9/20
        self.radius = radius
        self.slotwidth = slotwidth
        self.slotlength = slotlength
        self.xc = center
        self.xnw = [center[0] - 0.5*slotwidth,center[1] - (radius - slotlength)]
        self.xne = [center[0] + 0.5*slotwidth,center[1] - (radius - slotlength)]
        self.xsw = [center[0] - 0.5*slotwidth,center[1] - (radius)]
        self.xse = [center[0] + 0.5*slotwidth,center[1] - (radius)]
    def uOfXT(self,x,t):
        from math import sqrt
        dist = lambda u,v: sqrt( (u[0] - v[0])**2 + (u[1] - v[1])**2)
        phic = dist(self.xc,x) - self.radius
        phine = -dist(self.xne,x)
        phinw = -dist(self.xnw,x)
        phise = dist(self.xse,x)
        phisw = dist(self.xsw,x)
        phin = self.xnw[1] - x[1]
        phis = -(self.xsw[1] - x[1])
        phie = self.xne[0] - x[0]
        phiw = -(self.xnw[0] - x[0])
        if x[1] >= self.xnw[1]:
            if x[0] < self.xnw[0]:
                phi = max(phic,phinw)
            else:
                if x[0] < self.xne[0]:
                    phi = max(phic,phin)
                else:
                    phi = max(phic,phine)
        elif x[1] >= self.xsw[1]:
            if x[0] < self.xnw[0]:
                phi = max(phic,phiw)
            else:
                if x[0] < self.xne[0]:
                    phi = min([phin,phie,phiw])
                else:
                    phi = max(phic,phie)
        else:
            if x[0] < self.xsw[0]:
                phi = phic
            else:
                if x[0] < self.xse[0]:
                    phi = min(phisw,phise)
                else:
                    phi = phic
        return -phi

if pseudo2D:        
    analyticalSolution = {0:OscillatingVortex3Dcylinder(L)}
else:
    analyticalSolution = {0:OscillatingVortex3D(L)}

class UnitSquareVortex(TransportCoefficients.TC_base):
    from proteus.ctransportCoefficients import unitSquareVortexEvaluate
    from proteus.ctransportCoefficients import unitSquareVortexLevelSetEvaluate
    def __init__(self,useHJ=False,epsFact=1.5,checkMass=False):
        self.epsFact=epsFact
        self.useHJ = useHJ
        mass={0:{0:'linear'}}
        advection={0:{0:'linear'}}
        diffusion={}
        potential={}
        reaction={}
        if self.useHJ:
            hamiltonian={0:{0:'linear'}}
        else:
            hamiltonian={}
        TransportCoefficients.TC_base.__init__(self,
                                             1,
                                             mass,
                                             advection,
                                             diffusion,
                                             potential,
                                             reaction,
                                             hamiltonian)
        self.checkMass=checkMass
        self.useMetrics = 0.0
	self.sc_uref=1.0
	self.sc_beta=1.0
    def attachModels(self,modelList):
        import math
        self.model = modelList[0]
	self.u_old_dof = numpy.copy(self.model.u[0].dof)
        self.q_v = numpy.copy(self.model.q['x'])
        self.q_v *= -1.0
        ox = 0.5*L[0]
        oy = 0.5*L[1]
        self.ebqe_v = numpy.zeros(self.model.ebqe[('dH',0,0)].shape,'d')
        for eN in range(self.q_v.shape[0]):
            for k in range(self.q_v.shape[1]):
                self.q_v[eN,k,0] = 2.0*math.pi*(oy - self.model.q['x'][eN,k,1])
                self.q_v[eN,k,1] = -2.0*math.pi*(ox - self.model.q['x'][eN,k,0])
        for ebNE in range(self.ebqe_v.shape[0]):
            for kb in range(self.ebqe_v.shape[1]):
                self.ebqe_v[ebNE,kb,0] = 2.0*math.pi*(oy - self.model.ebqe['x'][ebNE,kb,1])
                self.ebqe_v[ebNE,kb,1] = -2.0*math.pi*(ox - self.model.ebqe['x'][ebNE,kb,0])
        self.model.q[('velocity',0)]=self.q_v
        self.model.ebqe[('velocity',0)]=self.ebqe_v
    def preStep(self,t,firstStep=False):
        copyInstructions = {}
        return copyInstructions
    def postStep(self,t,firstStep=False):
        copyInstructions = {}
        return copyInstructions
    def evaluate(self,t,c):
        pass
coefficients = UnitSquareVortex(useHJ=True,epsFact=epsFactHeaviside,checkMass=checkMass)

coefficients.variableNames=['u']

#now define the Dirichlet boundary conditions

def getDBC(x,flag):
    pass

dirichletConditions = {0:getDBC}

initialConditions  = {0:analyticalSolution[0]}

fluxBoundaryConditions = {0:'outFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

## @}
