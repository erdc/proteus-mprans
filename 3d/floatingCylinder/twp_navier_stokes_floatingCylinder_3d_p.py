"""
Incompressible Navier-Stokes flow around a square obstacle in 3D.
"""
from proteus import *
from proteus.default_p import *
import sys
from math import *
from floatingCylinder import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
from proteus.mprans import RANS2P

LevelModelType = RANS2P.LevelModel
coefficients = RANS2P.Coefficients(epsFact=epsFact_viscosity,
                                   sigma=sigma_01,
                                   rho_0=rho_0,
                                   nu_0=nu_0,
                                   rho_1=rho_1,
                                   nu_1=nu_1,
                                   g=g,
                                   nd=nd,
                                   LS_model=1,
                                   epsFact_density=epsFact_density,
                                   stokes=useStokes,
                                   useRBLES=0.0,
                                   useMetrics=0.0,
                                   movingDomain=movingDomain,
                                   useConstantH=useConstantH)

coefficients.waterLevel=waterLevel
import ode
class RigidCylinder(AuxiliaryVariables.AV_base):
    def __init__(self,rho=1.0,center=(0,0,0),radius=1.0,length=1.0):
        import pdb
        #pdb.set_trace()
        self.mass = length*0.5*rho_0*pi*radius**2 +length*0.5*rho_1*pi*radius**2
        self.world = ode.World()
        self.world.setGravity( tuple(g) )
        # Create a body inside the world
        self.body = ode.Body(self.world)
        self.M = ode.Mass()
        #pdb.set_trace()
        self.M.setCylinder(rho,direction=3,r=radius,h=length)
        self.body.setMass(self.M)
        self.body.setPosition(center)
        self.last_position=center
        self.position=center
        self.last_velocity=(0.0,0.0,0.0)
        self.velocity=(0.0,0.0,0.0)
        self.h=(0.0,0.0,0.0)
    def attachModel(self,model,ar):
        import copy
        self.model=model
        self.ar=ar
        self.writer = Archiver.XdmfWriter()
        self.nd = model.levelModelList[-1].nSpace_global
        m = self.model.levelModelList[-1]
        flagMax = max(m.mesh.elementBoundaryMaterialTypes)
        flagMin = min(m.mesh.elementBoundaryMaterialTypes)
        assert(flagMin == 0)
        assert(flagMax >= 0)
        self.nForces=flagMax+1
        self.levelFlist=[]
        for m in self.model.levelModelList:
            if self.nd == 2:
                F = numpy.zeros((self.nForces,2),'d')
            elif self.nd == 3:
                F = numpy.zeros((self.nForces,3),'d')
            else:
                logEvent("Can't use stress computation for nd = "+`self.nd`)
                F=None
            self.levelFlist.append(F)
        self.historyF=[]
        self.historyF.append(copy.deepcopy(self.levelFlist))
        return self
    def get_u(self):
        #print "obastacle-u",self.last_velocity[0]
        return self.last_velocity[0]
    def get_v(self):
        #print "obastacle-v",self.last_velocity[1]
        return self.last_velocity[1]
    def get_w(self):
        #print "obastacle-v",self.last_velocity[2]
        return self.last_velocity[2]
    def calculate(self):
        import pdb
        for m,F in zip(self.model.levelModelList,self.levelFlist):
            F.flat[:]=0.0
            assert(self.nd ==3)
            print "----------------need to add force calculation-------------------------------"
            # cfemIntegrals.calculateExteriorElementBoundaryStress3D(m.mesh.elementBoundaryMaterialTypes,
            #                                                        m.mesh.exteriorElementBoundariesArray,
            #                                                        m.mesh.elementBoundaryElementsArray,
            #                                                        m.mesh.elementBoundaryLocalElementBoundariesArray,
            #                                                        m.ebqe[('u',0)],#pressure
            #                                                        m.ebqe[('velocity',1)],#mom_flux_vec_u #cek todo, need to add real momentum flux
            #                                                        m.ebqe[('velocity',2)],#mom_flux_vec_v
            #                                                        m.ebqe[('velocity',3)],#mom_flux_vec_w
            #                                                        m.ebqe[('dS_u',0)],#dS
            #                                                        m.ebqe[('n')],
            #                                                        F)
            logEvent("Force")
            logEvent(`F`)
            Ftot=F[0,:]
            for ib in range(1,self.nForces):
                Ftot+=F[ib,:]
            logEvent("Total force on all boundaries")
            logEvent(`Ftot`)
        logEvent("x Force " +`self.model.stepController.t_model`+" "+`F[-1,0]`)
        logEvent("y Force " +`self.model.stepController.t_model`+" "+`F[-1,1]`)
        logEvent("z Force " +`self.model.stepController.t_model`+" "+`F[-1,2]`)
        #assume moving in the x direction
        self.body.addForce((0.0,0.0,F[-1,1]))#constrain to vertical motion initially
        self.world.step(self.model.stepController.dt_model)
        #f = m a = m (v_new - v_old)/dt
        #f dt/m = v_new - v_old
        #v_new = v_old + f dt/m 
        print "net acceleration====================",F[-1,2]/self.mass+g[2]
        self.velocity = (0.0,
                         0.0,
                         self.last_velocity[2]+F[-1,2]*self.model.stepController.dt_model/self.mass+g[2]*self.model.stepController.dt_model)
        self.position = (self.last_position[0]+self.velocity[0]*self.model.stepController.dt_model,
                         self.last_position[1]+self.velocity[1]*self.model.stepController.dt_model,
                         self.last_position[2]+self.velocity[2]*self.model.stepController.dt_model)
        self.h = (self.velocity[0]*self.model.stepController.dt_model,
                  self.velocity[1]*self.model.stepController.dt_model,
                  self.velocity[2]*self.model.stepController.dt_model)
        #x,y,z = self.body.getPosition()
        #u,v,w = self.body.getLinearVel()
        #self.position=(x,y,z)
        #self.velocity=(u,v,w)
        #self.h = (self.velocity[0]*self.model.stepController.dt_model,
        #          self.velocity[1]*self.model.stepController.dt_model,
        #          self.velocity[2]*self.model.stepController.dt_model)
        print "%1.2fsec: pos=(%6.3f, %6.3f, %6.3f)  vel=(%6.3f, %6.3f, %6.3f)" % \
            (self.model.stepController.t_model, 
             self.position[0], self.position[1], self.position[2], 
             self.velocity[0],self.velocity[1],self.velocity[2])
        print "%1.2fsec: last_pos=(%6.3f, %6.3f, %6.3f)  last_vel=(%6.3f, %6.3f, %6.3f)" % \
            (self.model.stepController.t_model, 
             self.last_position[0], self.last_position[1], self.last_position[2], 
             self.last_velocity[0],self.last_velocity[1],self.last_velocity[2])
        self.h = (self.velocity[0]*self.model.stepController.dt_model,
                  self.velocity[1]*self.model.stepController.dt_model,
                  self.velocity[2]*self.model.stepController.dt_model)
        self.last_velocity=self.velocity
        self.last_position=self.position
        print "dt model in object ",self.model.stepController.dt_model

rc = RigidCylinder(rho=0.5*rho_0,center=cylinder_center+(0.0,),radius=cylinder_radius,length=domain.L[1])

def velRamp(t):
    if t < residence_time:
        return 1.0-exp(-25.0*t/residence_time)
    else:
        return 1.0

#hack
#open on top
#closed on sides 
#no slip on obstacle
def getDBC_p(x,flag):
    if altBC:
        if flag in [boundaryTags['downstream'],boundaryTags['top']]:
            return lambda x,t: -coefficients.g[2]*(rho_0*(inflow_height - x[2])
                                                   -(rho_0-rho_1)*smoothedHeaviside_integral(epsFact_density*he,inflow_height-waterLevel)
                                                   +(rho_0-rho_1)*smoothedHeaviside_integral(epsFact_density*he,x[2]-waterLevel))
    else:
        if flag == boundaryTags['downstream']:
            return lambda x,t: -coefficients.g[2]*(rho_0*(inflow_height - x[2])
                                                   -(rho_0-rho_1)*smoothedHeaviside_integral(epsFact_density*he,inflow_height-waterLevel)
                                                   +(rho_0-rho_1)*smoothedHeaviside_integral(epsFact_density*he,x[2]-waterLevel))

def getDBC_u(x,flag):
    if altBC:
        if flag in [boundaryTags['upstream'],boundaryTags['bottom'],boundaryTags['top']]:
            return lambda x,t: Um
        elif flag in [boundaryTags['obstacle']]:
            return lambda x,t: 0.0
    else:
        if flag == boundaryTags['upstream']:
            return lambda x,t: Um

def getDBC_v(x,flag):
    if altBC:
        if flag in [boundaryTags['upstream'],boundaryTags['downstream'],boundaryTags['bottom'],boundaryTags['top']]:
            return lambda x,t: 0.0
        elif flag in [boundaryTags['obstacle']]:
            return lambda x,t: 0.0
    else:
        if flag == boundaryTags['upstream']:
            return lambda x,t: 0.0

def getDBC_w(x,flag):
    if altBC:
        if flag in [boundaryTags['upstream'],boundaryTags['downstream'],boundaryTags['bottom']]:
            return lambda x,t: 0.0
        elif flag in [boundaryTags['obstacle']]:
            return lambda x,t: 0.0
    else:
        if flag == boundaryTags['upstream']:
            return lambda x,t: 0.0

dirichletConditions = {0:getDBC_p,
                       1:getDBC_u,
                       2:getDBC_v,
                       3:getDBC_w}

def getAFBC_p(x,flag):
    if altBC:
        if flag in [boundaryTags['obstacle'],boundaryTags['bottom'],boundaryTags['front'],boundaryTags['back']]:
            return lambda x,t: 0.0
        if flag == boundaryTags['upstream']:
            return lambda x,t: -Um
    else:
        if flag in [boundaryTags['obstacle'],boundaryTags['top'],
                    boundaryTags['bottom'],boundaryTags['front'],boundaryTags['back']]:
            return lambda x,t: 0.0
        if flag == boundaryTags['upstream']:
            return lambda x,t: -Um
    if flag == 0:
        return lambda x,t: 0.0

def getAFBC_u(x,flag):
    if altBC:
        if flag in [boundaryTags['front'],boundaryTags['back']]:
            return lambda x,t: 0.0
    else:
        if flag in [boundaryTags['obstacle'],boundaryTags['top'],
                    boundaryTags['bottom'],boundaryTags['front'],boundaryTags['back']]:
            return lambda x,t: 0.0
    if flag == 0:
        return lambda x,t: 0.0

def getAFBC_v(x,flag):
    if altBC:
        if flag in [boundaryTags['front'],boundaryTags['back']]:
            return lambda x,t: 0.0
    else:
        if flag in [boundaryTags['obstacle'],boundaryTags['top'],
                    boundaryTags['bottom'],boundaryTags['front'],boundaryTags['back']]:
            return lambda x,t: 0.0
    if flag == 0:
        return lambda x,t: 0.0

def getAFBC_w(x,flag):
    if altBC:
        if flag in [boundaryTags['front'],boundaryTags['back']]:
            return lambda x,t: 0.0
    else:
        if flag in [boundaryTags['obstacle'],boundaryTags['top'],
                    boundaryTags['bottom'],boundaryTags['front'],boundaryTags['back']]:
            return lambda x,t: 0.0
    if flag == 0:
        return lambda x,t: 0.0

def getDFBC_u(x,flag):
    if altBC:
        #downstream is outflow
        if flag in [boundaryTags['front'],boundaryTags['back'],boundaryTags['downstream']]:
            return lambda x,t: 0.0
    else:
        if flag in [boundaryTags['obstacle'],boundaryTags['top'],boundaryTags['downstream'],
                    boundaryTags['bottom'],boundaryTags['front'],boundaryTags['back']]:
            return lambda x,t: 0.0
    if flag == 0:
        return lambda x,t: 0.0

def getDFBC_v(x,flag):
    if altBC:
        if flag in [boundaryTags['front'],boundaryTags['back']]:
            return lambda x,t: 0.0
    else:
        if flag in [boundaryTags['obstacle'],boundaryTags['top'],boundaryTags['downstream'],
                    boundaryTags['bottom'],boundaryTags['front'],boundaryTags['back']]:
            return lambda x,t: 0.0
    if flag == 0:
        return lambda x,t: 0.0

def getDFBC_w(x,flag):
    if altBC:
        #top is outflow
        if flag in [boundaryTags['top'],boundaryTags['front'],boundaryTags['back']]:
            return lambda x,t: 0.0
    else:
        if flag in [boundaryTags['obstacle'],boundaryTags['top'],boundaryTags['downstream'],
                    boundaryTags['bottom'],boundaryTags['front'],boundaryTags['back']]:
            return lambda x,t: 0.0
    if flag == 0:
        return lambda x,t: 0.0

fluxBoundaryConditions = {0:'mixedFlow',
                          1:'mixedFlow',
                          2:'mixedFlow',
                          3:'mixedFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p,
                                    1:getAFBC_u,
                                    2:getAFBC_v,
                                    3:getAFBC_w}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u},
                                   2:{2:getDFBC_v},
                                   3:{3:getDFBC_w}}

class Steady_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return -coefficients.g[2]*(rho_0*(inflow_height - x[2])
                                   -(rho_0-rho_1)*smoothedHeaviside_integral(epsFact_density*he,inflow_height-waterLevel)
                                   +(rho_0-rho_1)*smoothedHeaviside_integral(epsFact_density*he,x[2]-waterLevel))

class Steady_u:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return Um

class Steady_v:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

class Steady_w:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:Steady_p(),
                     1:Steady_u(),
                     2:Steady_v(),
                     3:Steady_w()}
