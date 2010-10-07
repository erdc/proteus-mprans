import proteus
from proteus.mprans.cVRANS2P import *

class Coefficients(proteus.TransportCoefficients.TC_base):
    """
    The coefficients for two incompresslble fluids governed by the Navier-Stokes equations and separated by a sharp interface represented by a level set function
    """
    from proteus.ctransportCoefficients import TwophaseNavierStokes_ST_LS_SO_2D_Evaluate
    from proteus.ctransportCoefficients import TwophaseNavierStokes_ST_LS_SO_3D_Evaluate
    from proteus.ctransportCoefficients import TwophaseNavierStokes_ST_LS_SO_2D_Evaluate_sd
    from proteus.ctransportCoefficients import TwophaseNavierStokes_ST_LS_SO_3D_Evaluate_sd
    from proteus.ctransportCoefficients import calculateWaveFunction3d_ref
    def __init__(self,
                 epsFact=1.5,
                 sigma=72.8,
                 rho_0=998.2,nu_0=1.004e-6,
                 rho_1=1.205,nu_1=1.500e-5,
                 g=[0.0,-9.8],
                 nd=2,
                 LS_model=None,
                 KN_model=None,
                 epsFact_density=None,
                 stokes=False,
                 sd=True,
                 movingDomain=False,
                 useRBLES=0.0,
                 meanGrainSize=0.01,
                 setParamsFunc=None,      #uses setParamsFunc if given
                 meanGrainSizeTypes=None, #otherwise can use element constant values
                 porosityTypes=None,
                 killNonlinearDrag=False,
                 waveFlag=None,
                 waveHeight=0.01,
                 waveCelerity=1.0,
                 waveFrequency=1.0,
                 waveNumber=2.0,
                 waterDepth=0.5,
                 Omega_s=[[0.45,0.55],[0.2,0.4],[0.0,1.0]],
                 epsFact_source=1.):
        self.useRBLES=useRBLES
        self.sd=sd
        if epsFact_density != None:
            self.epsFact_density = epsFact_density
        else:
            self.epsFact_density = epsFact
        self.stokes=stokes
        self.LS_model=LS_model
        self.KN_model=KN_model
        self.epsFact=epsFact
        self.eps=None
        self.sigma=sigma
        self.rho_0 = rho_0
        self.nu_0 = nu_0
        #cek for debugging using single phase test problems
        self.rho=rho_0
        self.nu=nu_0
        self.rho_1 = rho_1
        self.nu_1 = nu_1
        self.g = numpy.array(g)
        self.nd=nd
        #
        self.meanGrainSize     = meanGrainSize
        self.setParamsFunc=setParamsFunc
        self.meanGrainSizeTypes = meanGrainSizeTypes
        self.porosityTypes      = porosityTypes
        self.killNonlinearDrag  = int(killNonlinearDrag)
        self.waveFlag=waveFlag
        self.waveHeight=waveHeight
        self.waveCelerity=waveCelerity
        self.waveFrequency=waveFrequency
        self.waveNumber=waveNumber
        self.waterDepth=waterDepth
        self.Omega_s=Omega_s
        self.epsFact_source=epsFact_source
        self.linearDragFactor = 1.0; self.nonlinearDragFactor = 1.0
        if self.killNonlinearDrag:
            self.nonlinearDragFactor = 0.0
        
        mass={}
        advection={}
        diffusion={}
        potential={}
        reaction={}
        hamiltonian={}
        if nd==2:
            variableNames=['p','u','v']
            mass= {1:{1:'linear'},
                   2:{2:'linear'}}
            advection = {0:{0:'linear',
                            1:'linear',
                            2:'linear'},
                         1:{0:'nonlinear',
                            1:'nonlinear',
                            2:'nonlinear'},
                         2:{0:'nonlinear',
                            1:'nonlinear',
                            2:'nonlinear'}}
            diffusion  = {1:{1:{1:'constant'},2:{2:'constant'}},
                          2:{2:{2:'constant'},1:{1:'constant'}}}
            sdInfo  = {(1,1):(numpy.array([0,1,2],dtype='i'),
                             numpy.array([0,1],dtype='i')),
                       (1,2):(numpy.array([0,0,1],dtype='i'),
                              numpy.array([0],dtype='i')),
                       (2,2):(numpy.array([0,1,2],dtype='i'),
                              numpy.array([0,1],dtype='i')),
                       (2,1):(numpy.array([0,1,1],dtype='i'),
                              numpy.array([1],dtype='i'))}
            potential= {1:{1:'u'},
                        2:{2:'u'}}
            reaction = {0:{0:'constant'},
                        1:{1:'nonlinear',2:'nonlinear'},
                        2:{1:'nonlinear',2:'nonlinear'}}
            hamiltonian = {1:{0:'linear'},
                           2:{0:'linear'}}
            TC_base.__init__(self,
                             3,
                             mass,
                             advection,
                             diffusion,
                             potential,
                             reaction,
                             hamiltonian,
                             variableNames,
                             sparseDiffusionTensors=sdInfo,
                             useSparseDiffusion = sd,
                             movingDomain=movingDomain)
            self.vectorComponents=[1,2]
        elif nd==3:
            variableNames=['p','u','v','w']
            mass = {1:{1:'linear'},
                    2:{2:'linear'},
                    3:{3:'linear'}}
            advection = {0:{1:'linear',
                            2:'linear',
                            3:'linear'},
                         1:{0:'nonlinear',
                            1:'nonlinear',
                            2:'nonlinear',
                            3:'nonlinear'},
                         2:{0:'nonlinear',
                            1:'nonlinear',
                            2:'nonlinear',
                            3:'nonlinear'},
                         3:{0:'nonlinear',
                            1:'nonlinear',
                            2:'nonlinear',
                            3:'nonlinear'}}
            diffusion = {1:{1:{1:'constant'},2:{2:'constant'},3:{3:'constant'}},
                         2:{1:{1:'constant'},2:{2:'constant'},3:{3:'constant'}},
                         3:{1:{1:'constant'},2:{2:'constant'},3:{3:'constant'}}}
            sdInfo={}
            sdInfo  = {(1,1):(numpy.array([0,1,2,3],dtype='i'),numpy.array([0,1,2],dtype='i')),
                       (1,2):(numpy.array([0,0,1,1],dtype='i'),numpy.array([0],dtype='i')),
                       (1,3):(numpy.array([0,0,0,1],dtype='i'),numpy.array([0],dtype='i')),
                       (2,1):(numpy.array([0,1,1,1],dtype='i'),numpy.array([1],dtype='i')),
                       (2,2):(numpy.array([0,1,2,3],dtype='i'),numpy.array([0,1,2],dtype='i')),
                       (2,3):(numpy.array([0,0,0,1],dtype='i'),numpy.array([1],dtype='i')),
                       (3,1):(numpy.array([0,1,1,1],dtype='i'),numpy.array([2],dtype='i')),
                       (3,2):(numpy.array([0,0,1,1],dtype='i'),numpy.array([2],dtype='i')),
                       (3,3):(numpy.array([0,1,2,3],dtype='i'),numpy.array([0,1,2],dtype='i'))}
            potential= {1:{1:'u'},
                        2:{2:'u'},
                        3:{3:'u'}}
            reaction = {0:{0:'constant'},
                        1:{1:'nonlinear',2:'nonlinear',3:'nonlinear'},
                        2:{1:'nonlinear',2:'nonlinear',3:'nonlinear'},
                        3:{1:'nonlinear',2:'nonlinear',3:'nonlinear'}}
            hamiltonian = {1:{0:'linear'},
                           2:{0:'linear'},
                           3:{0:'linear'}}
            TC_base.__init__(self,
                             4,
                             mass,
                             advection,
                             diffusion,
                             potential,
                             reaction,
                             hamiltonian,
                             variableNames,
                             sparseDiffusionTensors=sdInfo,
                             useSparseDiffusion = sd,
                             movingDomain=movingDomain)
            self.vectorComponents=[1,2,3]
        
    def attachModels(self,modelList):
        #level set
        self.model = modelList[0]
        if self.LS_model != None:
            self.q_phi = modelList[self.LS_model].q[('u',0)]
            if modelList[self.LS_model].ebq.has_key(('u',0)):
                self.ebq_phi = modelList[self.LS_model].ebq[('u',0)]
            else:
                self.ebq_phi = None
            self.ebqe_phi = modelList[self.LS_model].ebqe[('u',0)]
            #normal
            self.q_n = modelList[self.LS_model].q[('grad(u)',0)]
            if modelList[self.LS_model].ebq.has_key(('grad(u)',0)):
                self.ebq_n = modelList[self.LS_model].ebq[('grad(u)',0)]
            else:
                self.ebq_n   = None
            self.ebqe_n    = modelList[self.LS_model].ebqe[('grad(u)',0)]
        #curvature
        if self.KN_model != None:
            self.q_kappa    = modelList[self.KN_model].q[('u',0)]
            self.ebqe_kappa = modelList[self.KN_model].ebqe[('u',0)]
            if modelList[self.KN_model].ebq.has_key(('u',0)):
                self.ebq_kappa = modelList[self.KN_model].ebq[('u',0)]
            else:
                self.ebq_kappa = None
    def initializeMesh(self,mesh):
        #cek we eventually need to use the local element diameter
        self.eps_density = self.epsFact_density*mesh.h
        self.eps_viscosity = self.epsFact*mesh.h
        self.mesh = mesh
        self.elementMaterialTypes = mesh.elementMaterialTypes
        self.eps_source=self.epsFact_source*mesh.h
    #initialize so it can run as single phase
    def initializeElementQuadrature(self,t,cq):
        if self.LS_model == None:
            self.q_phi = -numpy.ones(cq[('u',1)].shape,'d')
            self.q_n = -numpy.ones(cq[('velocity',0)].shape,'d')
        if self.KN_model == None:
            self.q_kappa = -numpy.zeros(cq[('u',1)].shape,'d')
        #VRANS
        self.q_porosity = numpy.ones(cq[('u',1)].shape,'d')
        self.q_meanGrain= numpy.ones(cq[('u',1)].shape,'d')            
        self.q_meanGrain.fill(self.meanGrainSize)
        if self.setParamsFunc != None:
            self.setParamsFunc(cq['x'],self.q_porosity,self.q_meanGrain)
        else:
            #TODO make loops faster
            if self.porosityTypes != None:
                for eN in range(self.q_porosity.shape[0]):
                    self.q_porosity[eN,:] = self.porosityTypes[self.elementMaterialTypes[eN]]
            if self.meanGrainSizeTypes != None:
                for eN in range(self.q_meanGrain.shape[0]):
                    self.q_meanGrain[eN,:] = self.meanGrainSizeTypes[self.elementMaterialTypes[eN]]
        #
    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        if self.LS_model == None:
            self.ebq_phi = -numpy.ones(cebq[('u',1)].shape,'d')
            self.ebq_n = -numpy.ones(cebq[('velocity',0)].shape,'d')
        if self.KN_model == None:
            self.ebq_kappa = -numpy.zeros(cebq[('u',1)].shape,'d')
        #VRANS
        self.ebq_porosity = numpy.ones(cebq[('u',1)].shape,'d')
        self.ebq_meanGrain= numpy.ones(cebq[('u',1)].shape,'d')
        self.ebq_meanGrain.fill(self.meanGrainSize)
        if self.setParamsFunc != None:
            self.setParamsFunc(cebq['x'],self.ebq_porosity,self.ebq_meanGrain)
        #TODO which mean to use or leave discontinuous
        #TODO make loops faster
        if self.porosityTypes != None:
            for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
                ebN = self.mesh.interiorElementBoundariesArray[ebNI]
                eN_left  = self.mesh.elementBoundaryElementsArray[ebN,0]
                eN_right = self.mesh.elementBoundaryElementsArray[ebN,1]
                ebN_element_left = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
                ebN_element_right = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,1]
                avg = 0.5*(self.porosityTypes[self.elementMaterialTypes[eN_left]]+
                           self.porosityTypes[self.elementMaterialTypes[eN_right]])
                self.ebq_porosity[eN_left,ebN_element_left,:]  = self.porosityTypes[self.elementMaterialTypes[eN_left]]
                self.ebq_porosity[eN_right,ebN_element_right,:]= self.porosityTypes[self.elementMaterialTypes[eN_right]]
            for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                eN  = self.mesh.elementBoundaryElementsArray[ebN,0]
                ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
                self.ebq_porosity[eN,ebN_element,:] = self.porosityTypes[self.elementMaterialTypes[eN]]
        if self.meanGrainSizeTypes != None:
            for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
                ebN = self.mesh.interiorElementBoundariesArray[ebNI]
                eN_left  = self.mesh.elementBoundaryElementsArray[ebN,0]
                eN_right = self.mesh.elementBoundaryElementsArray[ebN,1]
            ebN_element_left = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
            ebN_element_right = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,1]
            avg = 0.5*(self.meanGrainSizeTypes[self.elementMaterialTypes[eN_left]]+
                       self.meanGrainSizeTypes[self.elementMaterialTypes[eN_right]])
            self.ebq_meanGrain[eN_left,ebN_element_left,:] = self.meanGrainSizeTypes[self.elementMaterialTypes[eN_left]]
            self.ebq_meanGrain[eN_right,ebN_element_right,:] = self.meanGrainSizeTypes[self.elementMaterialTypes[eN_right]]
            for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                eN  = self.mesh.elementBoundaryElementsArray[ebN,0]
                ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
                self.ebq_meanGrain[eN,ebN_element,:] = self.meanGrainSizeTypes[self.elementMaterialTypes[eN]]
         #
    def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
        if self.LS_model == None:
            self.ebqe_phi = -numpy.ones(cebqe[('u',1)].shape,'d')
            self.ebqe_n = -numpy.ones(cebqe[('velocity',0)].shape,'d')
        if self.KN_model == None:
            self.ebqe_kappa = -numpy.zeros(cebqe[('u',1)].shape,'d')
        #VRANS
        self.ebqe_porosity = numpy.ones(cebqe[('u',1)].shape,'d')
        self.ebqe_meanGrain = numpy.ones(cebqe[('u',1)].shape,'d')
        self.ebqe_meanGrain.fill(self.meanGrainSize)
        #TODO make loops faster
        if self.setParamsFunc != None:
            self.setParamsFunc(cebqe['x'],self.ebqe_porosity,self.ebqe_meanGrain)
        else:
            if self.porosityTypes != None:
                for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                    ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                    eN  = self.mesh.elementBoundaryElementsArray[ebN,0]
                    self.ebqe_porosity[ebNE,:] = self.porosityTypes[self.elementMaterialTypes[eN]]
            if self.meanGrainSizeTypes != None:
                for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                    ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                    eN  = self.mesh.elementBoundaryElementsArray[ebN,0]
                    self.ebqe_meanGrain[ebNE,:] = self.meanGrainSizeTypes[self.elementMaterialTypes[eN]]
        #
    def updateToMovingDomain(self,t,c):
        from proteus import cfemIntegrals
        assert(self.movingDomain)
        if self.movingDomain:
            cfemIntegrals.update_f_movingDomain_constantMass(c['xt'],c[('f',0)])
            cfemIntegrals.update_f_movingDomain(c['xt'],c[('m',1)],c[('f',1)])
            cfemIntegrals.update_df_movingDomain(c['xt'],c[('dm',1,1)],c[('df',1,1)])
            cfemIntegrals.update_f_movingDomain(c['xt'],c[('m',2)],c[('f',2)])
            cfemIntegrals.update_df_movingDomain(c['xt'],c[('dm',2,2)],c[('df',2,2)])
            if self.nd == 3:
                cfemIntegrals.update_f_movingDomain(c['xt'],c[('m',3)],c[('f',3)])
                cfemIntegrals.update_df_movingDomain(c['xt'],c[('dm',3,3)],c[('df',3,3)])
    def evaluateForcingTerms(self,t,c,mesh=None,mesh_trial_ref=None,mesh_l2g=None):
        if c.has_key('x') and len(c['x'].shape) == 3:
            if self.nd == 2:
                #mwf debug
                #import pdb
                #pdb.set_trace()
                c[('r',0)].fill(0.0)
                eps_source=self.eps_source
                if self.waveFlag == 1:#secondOrderStokes:
                    waveFunctions.secondOrderStokesWave(c[('r',0)].shape[0],
                                                        c[('r',0)].shape[1],
                                                        self.waveHeight,
                                                        self.waveCelerity,
                                                        self.waveFrequency,
                                                        self.waveNumber,
                                                        self.waterDepth,
                                                        self.Omega_s[0][0],
                                                        self.Omega_s[0][1],
                                                        self.Omega_s[1][0],
                                                        self.Omega_s[1][1],
                                                        eps_source,
                                                        c['x'],
                                                        c[('r',0)],
                                                        t)
                elif self.waveFlag == 2:#solitary wave
                    waveFunctions.solitaryWave(c[('r',0)].shape[0],
                                               c[('r',0)].shape[1],
                                               self.waveHeight,
                                               self.waveCelerity,
                                               self.waveFrequency,
                                               self.waterDepth,
                                               self.Omega_s[0][0],
                                               self.Omega_s[0][1],
                                               self.Omega_s[1][0],
                                               self.Omega_s[1][1],
                                               eps_source,
                                               c['x'],
                                               c[('r',0)],
                                               t)

                elif self.waveFlag == 0:
                    waveFunctions.monochromaticWave(c[('r',0)].shape[0],
                                                    c[('r',0)].shape[1],
                                                    self.waveHeight,
                                                    self.waveCelerity,
                                                    self.waveFrequency,
                                                    self.Omega_s[0][0],
                                                    self.Omega_s[0][1],
                                                    self.Omega_s[1][0],
                                                    self.Omega_s[1][1],
                                                    eps_source,
                                                    c['x'],
                                                    c[('r',0)],
                                                    t)

                #mwf debug
                if numpy.isnan(c[('r',0)].any()):
                    import pdb
                    pdb.set_trace()
            else:
                #mwf debug
                #import pdb
                #pdb.set_trace()
                c[('r',0)].fill(0.0)
                eps_source=self.eps_source
                if self.waveFlag == 1:#secondOrderStokes:
                    waveFunctions.secondOrderStokesWave3d(c[('r',0)].shape[0],
                                                          c[('r',0)].shape[1],
                                                          self.waveHeight,
                                                          self.waveCelerity,
                                                          self.waveFrequency,
                                                          self.waveNumber,
                                                          self.waterDepth,
                                                          self.Omega_s[0][0],
                                                          self.Omega_s[0][1],
                                                          self.Omega_s[1][0],
                                                          self.Omega_s[1][1],
                                                          self.Omega_s[2][0],
                                                          self.Omega_s[2][1],
                                                          eps_source,
                                                          c['x'],
                                                          c[('r',0)],
                                                          t)
                elif self.waveFlag == 2:#solitary wave
                    waveFunctions.solitaryWave3d(c[('r',0)].shape[0],
                                                 c[('r',0)].shape[1],
                                                 self.waveHeight,
                                                 self.waveCelerity,
                                                 self.waveFrequency,
                                                 self.waterDepth,
                                                 self.Omega_s[0][0],
                                                 self.Omega_s[0][1],
                                                 self.Omega_s[1][0],
                                                 self.Omega_s[1][1],
                                                 self.Omega_s[2][0],
                                                 self.Omega_s[2][1],
                                                 eps_source,
                                                 c['x'],
                                                 c[('r',0)],
                                                 t)

                elif self.waveFlag == 0:
                    waveFunctions.monochromaticWave3d(c[('r',0)].shape[0],
                                                      c[('r',0)].shape[1],
                                                      self.waveHeight,
                                                      self.waveCelerity,
                                                      self.waveFrequency,
                                                      self.Omega_s[0][0],
                                                      self.Omega_s[0][1],
                                                      self.Omega_s[1][0],
                                                      self.Omega_s[1][1],
                                                      self.Omega_s[2][0],
                                                      self.Omega_s[2][1],
                                                      eps_source,
                                                      c['x'],
                                                      c[('r',0)],
                                                      t)

        else:
            assert mesh != None
            assert mesh_trial_ref != None
            assert mesh_l2g != None
            self.calculateWaveFunction3d_ref(mesh_trial_ref,
                                             mesh.nodeArray,
                                             mesh_l2g,
                                             mesh.elementDiametersArray,
                                             numpy.array(self.Omega_s[0]),
                                             numpy.array(self.Omega_s[1]),
                                             numpy.array(self.Omega_s[2]),
                                             t,
                                             self.waveFlag,
                                             self.epsFact_source,
                                             self.waveHeight,
                                             self.waveCelerity,
                                             self.waveFrequency,
                                             self.waveNumber,
                                             self.waterDepth,
                                             c[('r',0)])

    def evaluate(self,t,c):
        import pdb
        phi = None; n = None; kappa = None; porosity = None; meanGrain = None
 
        if c[('u',0)].shape == self.q_phi.shape:
            phi = self.q_phi
            n   = self.q_n
            kappa = self.q_kappa
            porosity = self.q_porosity
            meanGrain= self.q_meanGrain
        elif c[('u',0)].shape == self.ebqe_phi.shape:
            phi = self.ebqe_phi
            n = self.ebqe_n
            kappa = self.ebqe_kappa
            porosity = self.ebqe_porosity
            meanGrain= self.ebqe_meanGrain
        else:
            phi = self.ebq_phi
            n = self.ebq_n
            kappa = self.ebq_kappa
            porosity = self.ebq_porosity
            meanGrain= self.ebq_meanGrain
        #
        #mwf debug
        if phi == None or n == None or kappa == None or porosity == None or meanGrain == None:
            pdb.set_trace()
        #pdb.set_trace()
        if self.nd==2:
            if self.sd:
                self.VolumeAveragedTwophaseNavierStokes_ST_LS_SO_2D_Evaluate_sd(self.killNonlinearDrag,
                                                                                self.eps_density,
                                                                                self.eps_viscosity,
                                                                                self.sigma,
                                                                                self.rho_0,
                                                                                self.nu_0,
                                                                                self.rho_1,
                                                                                self.nu_1,
                                                                                meanGrain,
                                                                                self.g,
                                                                                phi,
                                                                                n,
                                                                                kappa,
                                                                                c[('u',0)],
                                                                                c[('grad(u)',0)],
                                                                                c[('u',1)],
                                                                                c[('u',2)],
                                                                                porosity,
                                                                                c[('m',1)],
                                                                                c[('dm',1,1)],
                                                                                c[('m',2)],
                                                                                c[('dm',2,2)],
                                                                                c[('f',0)],
                                                                                c[('df',0,1)],
                                                                                c[('df',0,2)],
                                                                                c[('f',1)],
                                                                                c[('df',1,1)],
                                                                                c[('df',1,2)],
                                                                                c[('f',2)],
                                                                                c[('df',2,1)],
                                                                                c[('df',2,2)],
                                                                                c[('a',1,1)],
                                                                                c[('a',2,2)],
                                                                                c[('a',1,2)],
                                                                                c[('a',2,1)],
                                                                                c[('r',1)],
                                                                                c[('r',2)],
                                                                                c[('dr',1,1)],
                                                                                c[('dr',1,2)],
                                                                                c[('dr',2,1)],
                                                                                c[('dr',2,2)],
                                                                                c[('H',1)],
                                                                                c[('dH',1,0)],
                                                                                c[('H',2)],
                                                                                c[('dH',2,0)])
            else:
                self.VolumeAveragedTwophaseNavierStokes_ST_LS_SO_2D_Evaluate(self.killNonlinearDrag,
                                                                             self.eps_density,
                                                                             self.eps_viscosity,
                                                                             self.sigma,
                                                                             self.rho_0,
                                                                             self.nu_0,
                                                                             self.rho_1,
                                                                             self.nu_1,
                                                                             meanGrain,
                                                                             self.g,
                                                                             phi,
                                                                             n,
                                                                             kappa,
                                                                             c[('u',0)],
                                                                             c[('grad(u)',0)],
                                                                             c[('u',1)],
                                                                             c[('u',2)],
                                                                             porosity,
                                                                             c[('m',1)],
                                                                             c[('dm',1,1)],
                                                                             c[('m',2)],
                                                                             c[('dm',2,2)],
                                                                             c[('f',0)],
                                                                             c[('df',0,1)],
                                                                             c[('df',0,2)],
                                                                             c[('f',1)],
                                                                             c[('df',1,1)],
                                                                             c[('df',1,2)],
                                                                             c[('f',2)],
                                                                             c[('df',2,1)],
                                                                             c[('df',2,2)],
                                                                             c[('a',1,1)],
                                                                             c[('a',2,2)],
                                                                             c[('a',1,2)],
                                                                             c[('a',2,1)],
                                                                             c[('r',1)],
                                                                             c[('r',2)],
                                                                             c[('dr',1,1)],
                                                                             c[('dr',1,2)],
                                                                             c[('dr',2,1)],
                                                                             c[('dr',2,2)],
                                                                             c[('H',1)],
                                                                             c[('dH',1,0)],
                                                                             c[('H',2)],
                                                                             c[('dH',2,0)])
            
            #
            if self.stokes:
                c[('f',1)].flat[:]=0.0; c[('f',2)].flat[:]=0.0;
                c[('df',1,1)].flat[:]=0.0;c[('df',1,2)].flat[:]=0.0;
                c[('df',2,2)].flat[:]=0.0;c[('df',2,1)].flat[:]=0.0;
            #
        elif  self.nd==3:
            if self.sd:
                self.VolumeAveragedTwophaseNavierStokes_ST_LS_SO_3D_Evaluate_sd(self.killNonlinearDrag,
                                                                                self.eps_density,
                                                                                self.eps_viscosity,
                                                                                self.sigma,
                                                                                self.rho_0,
                                                                                self.nu_0,
                                                                                self.rho_1,
                                                                                self.nu_1,
                                                                                meanGrain,
                                                                                self.g,
                                                                                phi,
                                                                                n,
                                                                                kappa,
                                                                                c[('u',0)],
                                                                                c[('grad(u)',0)],
                                                                                c[('u',1)],
                                                                                c[('u',2)],
                                                                                c[('u',3)],
                                                                                porosity,
                                                                                c[('m',1)],
                                                                                c[('dm',1,1)],
                                                                                c[('m',2)],
                                                                                c[('dm',2,2)],
                                                                                c[('m',3)],
                                                                                c[('dm',3,3)],
                                                                                c[('f',0)],
                                                                                c[('df',0,1)],
                                                                                c[('df',0,2)],
                                                                                c[('df',0,3)],
                                                                                c[('f',1)],
                                                                                c[('df',1,1)],
                                                                                c[('df',1,2)],
                                                                                c[('df',1,3)],
                                                                                c[('f',2)],
                                                                                c[('df',2,1)],
                                                                                c[('df',2,2)],
                                                                                c[('df',2,3)],
                                                                                c[('f',3)],
                                                                                c[('df',3,1)],
                                                                                c[('df',3,2)],
                                                                                c[('df',3,3)],
                                                                                c[('a',1,1)],
                                                                                c[('a',2,2)],
                                                                                c[('a',3,3)],
                                                                                c[('a',1,2)],
                                                                                c[('a',1,3)],
                                                                                c[('a',2,1)],
                                                                                c[('a',2,3)],
                                                                                c[('a',3,1)],
                                                                                c[('a',3,2)],
                                                                                c[('r',1)],
                                                                                c[('r',2)],
                                                                                c[('r',3)],
                                                                                c[('dr',1,1)],
                                                                                c[('dr',1,2)],
                                                                                c[('dr',1,3)],
                                                                                c[('dr',2,1)],
                                                                                c[('dr',2,2)],
                                                                                c[('dr',2,3)],
                                                                                c[('dr',3,1)],
                                                                                c[('dr',3,2)],
                                                                                c[('dr',3,3)],
                                                                                c[('H',1)],
                                                                                c[('dH',1,0)],
                                                                                c[('H',2)],
                                                                                c[('dH',2,0)],
                                                                                c[('H',3)],
                                                                                c[('dH',3,0)])
            else:
                self.VolumeAveragedTwophaseNavierStokes_ST_LS_SO_3D_Evaluate(self.killNonlinearDrag,
                                                                             self.eps_density,
                                                                             self.eps_viscosity,
                                                                             self.sigma,
                                                                             self.rho_0,
                                                                             self.nu_0,
                                                                             self.rho_1,
                                                                             self.nu_1,
                                                                             meanGrain,
                                                                             self.g,
                                                                             phi,
                                                                             n,
                                                                             kappa,
                                                                             c[('u',0)],
                                                                             c[('grad(u)',0)],
                                                                             c[('u',1)],
                                                                             c[('u',2)],
                                                                             c[('u',3)],
                                                                             porosity,
                                                                             c[('m',1)],
                                                                             c[('dm',1,1)],
                                                                             c[('m',2)],
                                                                             c[('dm',2,2)],
                                                                             c[('m',3)],
                                                                             c[('dm',3,3)],
                                                                             c[('f',0)],
                                                                             c[('df',0,1)],
                                                                             c[('df',0,2)],
                                                                             c[('df',0,3)],
                                                                             c[('f',1)],
                                                                             c[('df',1,1)],
                                                                             c[('df',1,2)],
                                                                             c[('df',1,3)],
                                                                             c[('f',2)],
                                                                             c[('df',2,1)],
                                                                             c[('df',2,2)],
                                                                             c[('df',2,3)],
                                                                             c[('f',3)],
                                                                             c[('df',3,1)],
                                                                             c[('df',3,2)],
                                                                             c[('df',3,3)],
                                                                             c[('a',1,1)],
                                                                             c[('a',2,2)],
                                                                             c[('a',3,3)],
                                                                             c[('a',1,2)],
                                                                             c[('a',1,3)],
                                                                             c[('a',2,1)],
                                                                             c[('a',2,3)],
                                                                             c[('a',3,1)],
                                                                             c[('a',3,2)],
                                                                             c[('r',1)],
                                                                             c[('r',2)],
                                                                             c[('r',3)],
                                                                             c[('dr',1,1)],
                                                                             c[('dr',1,2)],
                                                                             c[('dr',1,3)],
                                                                             c[('dr',2,1)],
                                                                             c[('dr',2,2)],
                                                                             c[('dr',2,3)],
                                                                             c[('dr',3,1)],
                                                                             c[('dr',3,2)],
                                                                             c[('dr',3,3)],
                                                                             c[('H',1)],
                                                                             c[('dH',1,0)],
                                                                             c[('H',2)],
                                                                             c[('dH',2,0)],
                                                                             c[('H',3)],
                                                                             c[('dH',3,0)])

        if c.has_key('x') and len(c['x'].shape) == 3:
            if self.nd == 2:
                #mwf debug
                #import pdb
                #pdb.set_trace()
                c[('r',0)].fill(0.0)
                eps_source=self.eps_source
                if self.waveFlag == 1:#secondOrderStokes:
                    waveFunctions.secondOrderStokesWave(c[('r',0)].shape[0],
                                                        c[('r',0)].shape[1],
                                                        self.waveHeight,
                                                        self.waveCelerity,
                                                        self.waveFrequency,
                                                        self.waveNumber,
                                                        self.waterDepth,
                                                        self.Omega_s[0][0],
                                                        self.Omega_s[0][1],
                                                        self.Omega_s[1][0],
                                                        self.Omega_s[1][1],
                                                        eps_source,
                                                        c['x'],
                                                        c[('r',0)],
                                                        t)
                elif self.waveFlag == 2:#solitary wave
                    waveFunctions.solitaryWave(c[('r',0)].shape[0],
                                               c[('r',0)].shape[1],
                                               self.waveHeight,
                                               self.waveCelerity,
                                               self.waveFrequency,
                                               self.waterDepth,
                                               self.Omega_s[0][0],
                                               self.Omega_s[0][1],
                                               self.Omega_s[1][0],
                                               self.Omega_s[1][1],
                                               eps_source,
                                               c['x'],
                                               c[('r',0)],
                                               t)

                elif self.waveFlag == 0:
                    waveFunctions.monochromaticWave(c[('r',0)].shape[0],
                                                    c[('r',0)].shape[1],
                                                    self.waveHeight,
                                                    self.waveCelerity,
                                                    self.waveFrequency,
                                                    self.Omega_s[0][0],
                                                    self.Omega_s[0][1],
                                                    self.Omega_s[1][0],
                                                    self.Omega_s[1][1],
                                                    eps_source,
                                                    c['x'],
                                                    c[('r',0)],
                                                    t)

                #mwf debug
                if numpy.isnan(c[('r',0)].any()):
                    import pdb
                    pdb.set_trace()
            else:
                #mwf debug
                #import pdb
                #pdb.set_trace()
                c[('r',0)].fill(0.0)
                eps_source=self.eps_source
                if self.waveFlag == 1:#secondOrderStokes:
                    waveFunctions.secondOrderStokesWave3d(c[('r',0)].shape[0],
                                                          c[('r',0)].shape[1],
                                                          self.waveHeight,
                                                          self.waveCelerity,
                                                          self.waveFrequency,
                                                          self.waveNumber,
                                                          self.waterDepth,
                                                          self.Omega_s[0][0],
                                                          self.Omega_s[0][1],
                                                          self.Omega_s[1][0],
                                                          self.Omega_s[1][1],
                                                          self.Omega_s[2][0],
                                                          self.Omega_s[2][1],
                                                          eps_source,
                                                          c['x'],
                                                          c[('r',0)],
                                                          t)
                elif self.waveFlag == 2:#solitary wave
                    waveFunctions.solitaryWave3d(c[('r',0)].shape[0],
                                                 c[('r',0)].shape[1],
                                                 self.waveHeight,
                                                 self.waveCelerity,
                                                 self.waveFrequency,
                                                 self.waterDepth,
                                                 self.Omega_s[0][0],
                                                 self.Omega_s[0][1],
                                                 self.Omega_s[1][0],
                                                 self.Omega_s[1][1],
                                                 self.Omega_s[2][0],
                                                 self.Omega_s[2][1],
                                                 eps_source,
                                                 c['x'],
                                                 c[('r',0)],
                                                 t)

                elif self.waveFlag == 0:
                    waveFunctions.monochromaticWave3d(c[('r',0)].shape[0],
                                                      c[('r',0)].shape[1],
                                                      self.waveHeight,
                                                      self.waveCelerity,
                                                      self.waveFrequency,
                                                      self.Omega_s[0][0],
                                                      self.Omega_s[0][1],
                                                      self.Omega_s[1][0],
                                                      self.Omega_s[1][1],
                                                      self.Omega_s[2][0],
                                                      self.Omega_s[2][1],
                                                      eps_source,
                                                      c['x'],
                                                      c[('r',0)],
                                                      t)

class LevelModel(proteus.Transport.OneLevelTransport):
    nCalls=0
    def __init__(self,
                 uDict,
                 phiDict,
                 testSpaceDict,
                 matType,
                 dofBoundaryConditionsDict,
                 dofBoundaryConditionsSetterDict,
                 coefficients,
                 elementQuadrature,
                 elementBoundaryQuadrature,
                 fluxBoundaryConditionsDict=None,
                 advectiveFluxBoundaryConditionsSetterDict=None,
                 diffusiveFluxBoundaryConditionsSetterDictDict=None,
                 stressTraceBoundaryConditionsSetterDictDict=None,
                 stabilization=None,
                 shockCapturing=None,
                 conservativeFluxDict=None,
                 numericalFluxType=None,
                 TimeIntegrationClass=None,
                 massLumping=False,
                 reactionLumping=False,
                 options=None,
                 name='VRANS2P',
                 reuse_trial_and_test_quadrature=True,
                 sd = True,
                 movingDomain=False):
        #
        #set the objects describing the method and boundary conditions
        #
        self.movingDomain=movingDomain
        self.tLast_mesh=None
        #
        #cek todo clean up these flags in the optimized version
        self.bcsTimeDependent=options.bcsTimeDependent
        self.bcsSet=False
        self.name=name
        self.sd=sd
        self.lowmem=True
        self.timeTerm=True#allow turning off  the  time derivative
        self.testIsTrial=True
        self.phiTrialIsTrial=True            
        self.u = uDict
        self.Hess=False
        if isinstance(self.u[0].femSpace,C0_AffineQuadraticOnSimplexWithNodalBasis):
            self.Hess=True
        self.ua = {}#analytical solutions
        self.phi  = phiDict
        self.dphi={}
        self.matType = matType
        #mwf try to reuse test and trial information across components if spaces are the same
        self.reuse_test_trial_quadrature = reuse_trial_and_test_quadrature#True#False
        if self.reuse_test_trial_quadrature:
            for ci in range(1,coefficients.nc):
                assert self.u[ci].femSpace.__class__.__name__ == self.u[0].femSpace.__class__.__name__, "to reuse_test_trial_quad all femSpaces must be the same!"
        ## Simplicial Mesh
        self.mesh = self.u[0].femSpace.mesh #assume the same mesh for  all components for now
        self.testSpace = testSpaceDict
        self.dirichletConditions = dofBoundaryConditionsDict
        self.dirichletNodeSetList=None #explicit Dirichlet  conditions for now, no Dirichlet BC constraints
        self.coefficients = coefficients
        self.coefficients.initializeMesh(self.mesh)
        self.nc = self.coefficients.nc
        self.stabilization = stabilization
        self.shockCapturing = shockCapturing
        self.conservativeFlux = conservativeFluxDict #no velocity post-processing for now
        self.fluxBoundaryConditions=fluxBoundaryConditionsDict
        self.advectiveFluxBoundaryConditionsSetterDict=advectiveFluxBoundaryConditionsSetterDict
        self.diffusiveFluxBoundaryConditionsSetterDictDict = diffusiveFluxBoundaryConditionsSetterDictDict
        #determine whether  the stabilization term is nonlinear
        self.stabilizationIsNonlinear = False
        #cek come back
	if self.stabilization != None:
	    for ci in range(self.nc):
		if coefficients.mass.has_key(ci):
		    for flag in coefficients.mass[ci].values():
			if flag == 'nonlinear':
			    self.stabilizationIsNonlinear=True
		if  coefficients.advection.has_key(ci):
		    for  flag  in coefficients.advection[ci].values():
			if flag == 'nonlinear':
			    self.stabilizationIsNonlinear=True
		if  coefficients.diffusion.has_key(ci):
		    for diffusionDict in coefficients.diffusion[ci].values():
			for  flag  in diffusionDict.values():
			    if flag != 'constant':
				self.stabilizationIsNonlinear=True
		if  coefficients.potential.has_key(ci):
 		    for flag in coefficients.potential[ci].values():
			if  flag == 'nonlinear':
			    self.stabilizationIsNonlinear=True
		if coefficients.reaction.has_key(ci):
		    for flag in coefficients.reaction[ci].values():
			if  flag == 'nonlinear':
			    self.stabilizationIsNonlinear=True
		if coefficients.hamiltonian.has_key(ci):
		    for flag in coefficients.hamiltonian[ci].values():
			if  flag == 'nonlinear':
			    self.stabilizationIsNonlinear=True
        #determine if we need element boundary storage
        self.elementBoundaryIntegrals = {}
        for ci  in range(self.nc):
            self.elementBoundaryIntegrals[ci] = ((self.conservativeFlux != None) or 
                                                 (numericalFluxType != None) or 
                                                 (self.fluxBoundaryConditions[ci] == 'outFlow') or
                                                 (self.fluxBoundaryConditions[ci] == 'mixedFlow') or
                                                 (self.fluxBoundaryConditions[ci] == 'setFlow'))
	#
        #calculate some dimensions
        #
        self.nSpace_global    = self.u[0].femSpace.nSpace_global #assume same space dim for all variables
        self.nDOF_trial_element     = [u_j.femSpace.max_nDOF_element for  u_j in self.u.values()]
        self.nDOF_phi_trial_element     = [phi_k.femSpace.max_nDOF_element for  phi_k in self.phi.values()]
        self.n_phi_ip_element = [phi_k.femSpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints for  phi_k in self.phi.values()]
        self.nDOF_test_element     = [femSpace.max_nDOF_element for femSpace in self.testSpace.values()]
        self.nFreeDOF_global  = [dc.nFreeDOF_global for dc in self.dirichletConditions.values()]
        self.nVDOF_element    = sum(self.nDOF_trial_element)
        self.nFreeVDOF_global = sum(self.nFreeDOF_global) 
        #
        NonlinearEquation.__init__(self,self.nFreeVDOF_global)
        #
        #build the quadrature point dictionaries from the input (this
        #is just for convenience so that the input doesn't have to be
        #complete)
        #
        elementQuadratureDict={}
        elemQuadIsDict = isinstance(elementQuadrature,dict)
        if elemQuadIsDict: #set terms manually
            for I in self.coefficients.elementIntegralKeys:
                if elementQuadrature.has_key(I):
                    elementQuadratureDict[I] = elementQuadrature[I]
                else:
                    elementQuadratureDict[I] = elementQuadrature['default']
        else:
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[I] = elementQuadrature
        if self.stabilization != None:
            for I in self.coefficients.elementIntegralKeys:
                if elemQuadIsDict:
                    if elementQuadrature.has_key(I):
                        elementQuadratureDict[('stab',)+I[1:]] = elementQuadrature[I]
                    else:
                        elementQuadratureDict[('stab',)+I[1:]] = elementQuadrature['default']
                else:
                    elementQuadratureDict[('stab',)+I[1:]] = elementQuadrature
        if self.shockCapturing != None:
            for ci in self.shockCapturing.components:
                if elemQuadIsDict:
                    if elementQuadrature.has_key(('numDiff',ci,ci)):
                        elementQuadratureDict[('numDiff',ci,ci)] = elementQuadrature[('numDiff',ci,ci)]
                    else:
                        elementQuadratureDict[('numDiff',ci,ci)] = elementQuadrature['default']
                else:
                    elementQuadratureDict[('numDiff',ci,ci)] = elementQuadrature
        if massLumping:
            for ci in self.coefficients.mass.keys():
                elementQuadratureDict[('m',ci)] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global,1)
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[('stab',)+I[1:]] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global,1)
        if reactionLumping:
            for ci in self.coefficients.mass.keys():
                elementQuadratureDict[('r',ci)] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global,1)
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[('stab',)+I[1:]] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global,1)
        elementBoundaryQuadratureDict={}
        if isinstance(elementBoundaryQuadrature,dict): #set terms manually
            for I in self.coefficients.elementBoundaryIntegralKeys:
                if elementBoundaryQuadrature.has_key(I):
                    elementBoundaryQuadratureDict[I] = elementBoundaryQuadrature[I]
                else:
                    elementBoundaryQuadratureDict[I] = elementBoundaryQuadrature['default']
        else:
            for I in self.coefficients.elementBoundaryIntegralKeys: 
                elementBoundaryQuadratureDict[I] = elementBoundaryQuadrature
        #
        # find the union of all element quadrature points and
        # build a quadrature rule for each integral that has a
        # weight at each point in the union
        #mwf include tag telling me which indices are which quadrature rule?
        (self.elementQuadraturePoints,self.elementQuadratureWeights,
         self.elementQuadratureRuleIndeces) = Quadrature.buildUnion(elementQuadratureDict)
        self.nQuadraturePoints_element = self.elementQuadraturePoints.shape[0]
        self.nQuadraturePoints_global = self.nQuadraturePoints_element*self.mesh.nElements_global
        #
        #Repeat the same thing for the element boundary quadrature
        #
        (self.elementBoundaryQuadraturePoints,
         self.elementBoundaryQuadratureWeights,
         self.elementBoundaryQuadratureRuleIndeces) = Quadrature.buildUnion(elementBoundaryQuadratureDict)
        self.nElementBoundaryQuadraturePoints_elementBoundary = self.elementBoundaryQuadraturePoints.shape[0]
        self.nElementBoundaryQuadraturePoints_global = (self.mesh.nElements_global*
                                                        self.mesh.nElementBoundaries_element*
                                                        self.nElementBoundaryQuadraturePoints_elementBoundary)

#        if isinstance(self.u[0].femSpace,C0_AffineLinearOnSimplexWithNodalBasis):
#            print self.nQuadraturePoints_element
#            if self.nSpace_global == 3:
#                assert(self.nQuadraturePoints_element == 5)
#            elif self.nSpace_global == 2:
#                assert(self.nQuadraturePoints_element == 6)
#            elif self.nSpace_global == 1:
#                assert(self.nQuadraturePoints_element == 3)
#
#            print self.nElementBoundaryQuadraturePoints_elementBoundary
#            if self.nSpace_global == 3:
#                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 4)
#            elif self.nSpace_global == 2:
#                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 4)
#            elif self.nSpace_global == 1:
#                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 1)
        #
        #simplified allocations for test==trial and also check if space is mixed or not
        #
        self.q={}
        self.ebq={}
        self.ebq_global={}
        self.ebqe={}
        self.phi_ip={}
        #mesh
        self.ebqe['x'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,3),'d')
        self.ebq_global[('totalFlux',0)] = numpy.zeros((self.mesh.nElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebq_global[('velocityAverage',0)] = numpy.zeros((self.mesh.nElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        self.q[('u',1)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('u',2)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('u',3)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m',1)] = self.q[('u',1)]
        self.q[('m',2)] = self.q[('u',2)]
        self.q[('m',3)] = self.q[('u',3)]
        self.q[('m_last',1)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m_last',2)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m_last',3)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m_tmp',1)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m_tmp',2)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m_tmp',3)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('f',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),'d')
        self.q[('velocity',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),'d')
        self.q[('cfl',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('numDiff',1,1)] =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('numDiff',2,2)] =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('numDiff',3,3)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.ebqe[('u',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('u',1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('u',2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('u',3)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('advectiveFlux_bc_flag',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('advectiveFlux_bc_flag',1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('advectiveFlux_bc_flag',2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('advectiveFlux_bc_flag',3)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('diffusiveFlux_bc_flag',1,1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('diffusiveFlux_bc_flag',2,2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('diffusiveFlux_bc_flag',3,3)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('advectiveFlux_bc',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('advectiveFlux_bc',1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('advectiveFlux_bc',2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('advectiveFlux_bc',3)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('diffusiveFlux_bc',1,1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('penalty')] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('diffusiveFlux_bc',2,2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('diffusiveFlux_bc',3,3)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('velocity',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        self.ebqe[('velocity',1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        self.ebqe[('velocity',2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        self.ebqe[('velocity',3)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        #VRANS start, defaults to RANS 
        self.q[('r',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        #VRANS end
        self.points_elementBoundaryQuadrature= set()
        self.scalars_elementBoundaryQuadrature= set([('u',ci) for ci in range(self.nc)])
        self.vectors_elementBoundaryQuadrature= set()
        self.tensors_elementBoundaryQuadrature= set()
        #
        #show quadrature
        #
        log("Dumping quadrature shapes for model %s" % self.name,level=9)
        log("Element quadrature array (q)", level=9)
        for (k,v) in self.q.iteritems(): log(str((k,v.shape)),level=9)
        log("Element boundary quadrature (ebq)",level=9) 
        for (k,v) in self.ebq.iteritems(): log(str((k,v.shape)),level=9)
        log("Global element boundary quadrature (ebq_global)",level=9)
        for (k,v) in self.ebq_global.iteritems(): log(str((k,v.shape)),level=9)
        log("Exterior element boundary quadrature (ebqe)",level=9)
        for (k,v) in self.ebqe.iteritems(): log(str((k,v.shape)),level=9)
        log("Interpolation points for nonlinear diffusion potential (phi_ip)",level=9)
        for (k,v) in self.phi_ip.iteritems(): log(str((k,v.shape)),level=9)
        #
        # allocate residual and Jacobian storage
        #
	self.inflowBoundaryBC = {}
	self.inflowBoundaryBC_values = {}
	self.inflowFlux = {}
 	for cj in range(self.nc):
 	    self.inflowBoundaryBC[cj] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,),'i')
 	    self.inflowBoundaryBC_values[cj] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nDOF_trial_element[cj]),'d')
 	    self.inflowFlux[cj] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.internalNodes = set(range(self.mesh.nNodes_global))
	#identify the internal nodes this is ought to be in mesh
        ##\todo move this to mesh
        for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
            ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
            eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
            ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
            for i in range(self.mesh.nNodes_element):
                if i != ebN_element:
                    I = self.mesh.elementNodesArray[eN_global,i]
                    self.internalNodes -= set([I])
        self.nNodes_internal = len(self.internalNodes)
        self.internalNodesArray=numpy.zeros((self.nNodes_internal,),'i')
        for nI,n in enumerate(self.internalNodes):
            self.internalNodesArray[nI]=n
        #
        del self.internalNodes
        self.internalNodes = None
        log("Updating local to global mappings",2)
        self.updateLocal2Global()
        log("Building time integration object",2)
        log(memory("inflowBC, internalNodes,updateLocal2Global","OneLevelTransport"),level=4)
        #mwf for interpolating subgrid error for gradients etc
        if self.stabilization and self.stabilization.usesGradientStabilization:
            self.timeIntegration = TimeIntegrationClass(self,integrateInterpolationPoints=True)
        else:
             self.timeIntegration = TimeIntegrationClass(self)
           
        if options != None:
            self.timeIntegration.setFromOptions(options)
        log(memory("TimeIntegration","OneLevelTransport"),level=4)
        log("Calculating numerical quadrature formulas",2)
        self.calculateQuadrature()
        #lay out components/equations contiguously for now
        self.offset = [0]
	for ci in range(1,self.nc):
	    self.offset += [self.offset[ci-1]+self.nFreeDOF_global[ci-1]]
        self.stride = [1 for ci in range(self.nc)]
        #use contiguous layout of components for parallel, requires weak DBC's
        comm = Comm.get()
        self.comm=comm
        if comm.size() > 1:
            assert numericalFluxType != None and numericalFluxType.useWeakDirichletConditions,"You must use a numerical flux to apply weak boundary conditions for parallel runs"
            self.offset = [0]
            for ci in range(1,self.nc):
                self.offset += [ci]
            self.stride = [self.nc for ci in range(self.nc)]
        #
        log(memory("stride+offset","OneLevelTransport"),level=4)
        if numericalFluxType != None:
            if options == None or options.periodicDirichletConditions == None:
                self.numericalFlux = numericalFluxType(self,
                                                       dofBoundaryConditionsSetterDict,
                                                       advectiveFluxBoundaryConditionsSetterDict,
                                                       diffusiveFluxBoundaryConditionsSetterDictDict)
            else:
                self.numericalFlux = numericalFluxType(self,
                                                       dofBoundaryConditionsSetterDict,
                                                       advectiveFluxBoundaryConditionsSetterDict,
                                                       diffusiveFluxBoundaryConditionsSetterDictDict,
                                                       options.periodicDirichletConditions)
        else:
            self.numericalFlux = None
        #set penalty terms
        #cek todo move into numerical flux initialization
        if self.ebq_global.has_key('penalty'):
            for ebN in range(self.mesh.nElementBoundaries_global):
                for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                    self.ebq_global['penalty'][ebN,k] = self.numericalFlux.penalty_constant/(self.mesh.elementBoundaryDiametersArray[ebN]**self.numericalFlux.penalty_power)
        #penalty term
        #cek move  to Numerical flux initialization
        if self.ebqe.has_key('penalty'):
            for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                    self.ebqe['penalty'][ebNE,k] = self.numericalFlux.penalty_constant/self.mesh.elementBoundaryDiametersArray[ebN]**self.numericalFlux.penalty_power
        log(memory("numericalFlux","OneLevelTransport"),level=4)
        self.elementEffectiveDiametersArray  = self.mesh.elementInnerDiametersArray
        #use post processing tools to get conservative fluxes, None by default
        from proteus import PostProcessingTools
        self.velocityPostProcessor = PostProcessingTools.VelocityPostProcessingChooser(self)  
        log(memory("velocity postprocessor","OneLevelTransport"),level=4)
        #helper for writing out data storage
        from proteus import Archiver
        self.elementQuadratureDictionaryWriter = Archiver.XdmfWriter()
        self.elementBoundaryQuadratureDictionaryWriter = Archiver.XdmfWriter()
        self.exteriorElementBoundaryQuadratureDictionaryWriter = Archiver.XdmfWriter()
        for ci,fbcObject  in self.fluxBoundaryConditionsObjectsDict.iteritems():
            self.ebqe[('advectiveFlux_bc_flag',ci)] = numpy.zeros(self.ebqe[('advectiveFlux_bc',ci)].shape,'i')
            for t,g in fbcObject.advectiveFluxBoundaryConditionsDict.iteritems():
                if self.coefficients.advection.has_key(ci):
                    self.ebqe[('advectiveFlux_bc',ci)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
                    self.ebqe[('advectiveFlux_bc_flag',ci)][t[0],t[1]] = 1
            for ck,diffusiveFluxBoundaryConditionsDict in fbcObject.diffusiveFluxBoundaryConditionsDictDict.iteritems():
                self.ebqe[('diffusiveFlux_bc_flag',ck,ci)] = numpy.zeros(self.ebqe[('diffusiveFlux_bc',ck,ci)].shape,'i')
                for t,g in diffusiveFluxBoundaryConditionsDict.iteritems():
                    self.ebqe[('diffusiveFlux_bc',ck,ci)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
                    self.ebqe[('diffusiveFlux_bc_flag',ck,ci)][t[0],t[1]] = 1
        self.numericalFlux.setDirichletValues(self.ebqe)

        #cek/ido todo replace python loops in modules with optimized code if possible/necessary
        self.forceStrongConditions=False
        self.dirichletConditionsForceDOF = {}
        if self.forceStrongConditions:
            for cj in range(self.nc):
                self.dirichletConditionsForceDOF[cj] = DOFBoundaryConditions(self.u[cj].femSpace,dofBoundaryConditionsSetterDict[cj],weakDirichletConditions=False)
        compKernelFlag = 0
        self.rans2p = cVRANS2P_base(self.nSpace_global,
                                   self.nQuadraturePoints_element,
                                   self.u[0].femSpace.elementMaps.localFunctionSpace.dim,
                                   self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.dim,
                                   self.testSpace[0].referenceFiniteElement.localFunctionSpace.dim,
                                   self.nElementBoundaryQuadraturePoints_elementBoundary,
                                   compKernelFlag)

    def getResidual(self,u,r):
        """
        Calculate the element residuals and add in to the global residual
        """
        #Load the unknowns into the finite element dof
        self.timeIntegration.calculateCoefs()
        self.timeIntegration.calculateU(u)
        self.setUnknowns(self.timeIntegration.u)
        #cek todo put in logic to skip if BC's don't depend on t or u
        #hack
        if self.bcsTimeDependent or not self.bcsSet:
            self.bcsSet=True
            #Dirichlet boundary conditions
            self.numericalFlux.setDirichletValues(self.ebqe)
            #Flux boundary conditions
            for ci,fbcObject  in self.fluxBoundaryConditionsObjectsDict.iteritems():
                for t,g in fbcObject.advectiveFluxBoundaryConditionsDict.iteritems():
                    if self.coefficients.advection.has_key(ci):
                        self.ebqe[('advectiveFlux_bc',ci)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
                        self.ebqe[('advectiveFlux_bc_flag',ci)][t[0],t[1]] = 1
                for ck,diffusiveFluxBoundaryConditionsDict in fbcObject.diffusiveFluxBoundaryConditionsDictDict.iteritems():
                    for t,g in diffusiveFluxBoundaryConditionsDict.iteritems():
                        self.ebqe[('diffusiveFlux_bc',ck,ci)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
                        self.ebqe[('diffusiveFlux_bc_flag',ck,ci)][t[0],t[1]] = 1
        r.fill(0.0)
        self.Ct_sge = 4.0
        self.Cd_sge = 144.0
 
        #TODO how to request problem specific evaluations from coefficient class
        if 'evaluateForcingTerms' in dir(self.coefficients):
            self.coefficients.evaluateForcingTerms(self.timeIntegration.t,self.q,self.mesh,
                                                   self.u[0].femSpace.elementMaps.psi,self.mesh.elementNodesArray)
        if self.forceStrongConditions:
            for cj in range(len(self.dirichletConditionsForceDOF)):
                for dofN,g in self.dirichletConditionsForceDOF[cj].DOFBoundaryConditionsDict.iteritems():
                    self.u[cj].dof[dofN] = g(self.dirichletConditionsForceDOF[cj].DOFBoundaryPointDict[dofN],self.timeIntegration.t)
 
        self.rans2p.calculateResidual(#element
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[1].femSpace.psi,
            self.u[1].femSpace.grad_psi,
            self.u[1].femSpace.psi,
            self.u[1].femSpace.grad_psi,
            #element boundary
            self.u[0].femSpace.elementMaps.psi_trace,
            self.u[0].femSpace.elementMaps.grad_psi_trace,
            self.elementBoundaryQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[1].femSpace.psi_trace,
            self.u[1].femSpace.grad_psi_trace,
            self.u[1].femSpace.psi_trace,
            self.u[1].femSpace.grad_psi_trace,
            self.u[0].femSpace.elementMaps.boundaryNormals,
            self.u[0].femSpace.elementMaps.boundaryJacobians,
            #physics
            self.mesh.elementDiametersArray,
            self.stabilization.hFactor,
            self.mesh.nElements_global,
            self.coefficients.useRBLES,
            self.timeIntegration.alpha_bdf,
            self.coefficients.epsFact_density,
            self.coefficients.epsFact,
            self.coefficients.sigma,
            self.coefficients.rho_0,
            self.coefficients.nu_0,
            self.coefficients.rho_1,
            self.coefficients.nu_1,
            self.Ct_sge,
            self.Cd_sge,
            self.shockCapturing.shockCapturingFactor,
            #VRANS start
            self.coefficients.linearDragFactor,
            self.coefficients.nonlinearDragFactor,
            self.coefficients.q_porosity,
            self.coefficients.q_meanGrain,
            self.q[('r',0)],
            #VRANS end
            self.u[0].femSpace.dofMap.l2g,
            self.u[1].femSpace.dofMap.l2g,
            self.u[0].dof,
            self.u[1].dof,
            self.u[2].dof,
            self.u[3].dof,
            self.coefficients.g,
            self.coefficients.q_phi,
            self.coefficients.q_n,
            self.coefficients.q_kappa,
            self.timeIntegration.m_tmp[1],
            self.timeIntegration.m_tmp[2],
            self.timeIntegration.m_tmp[3],
            self.q[('f',0)],
            self.timeIntegration.beta_bdf[1],
            self.timeIntegration.beta_bdf[2],
            self.timeIntegration.beta_bdf[3],
            self.stabilization.v_last,
            self.q[('cfl',0)],
            self.q[('numDiff',1,1)], 
            self.q[('numDiff',2,2)], 
            self.q[('numDiff',3,3)],
            self.shockCapturing.numDiff_last[1],
            self.shockCapturing.numDiff_last[2],
            self.shockCapturing.numDiff_last[3],
            self.coefficients.sdInfo[(1,1)][0],self.coefficients.sdInfo[(1,1)][1],
            self.coefficients.sdInfo[(1,2)][0],self.coefficients.sdInfo[(1,2)][1],
            self.coefficients.sdInfo[(1,3)][0],self.coefficients.sdInfo[(1,3)][1],
            self.coefficients.sdInfo[(2,2)][0],self.coefficients.sdInfo[(2,2)][1],
            self.coefficients.sdInfo[(2,1)][0],self.coefficients.sdInfo[(2,1)][1],
            self.coefficients.sdInfo[(2,3)][0],self.coefficients.sdInfo[(2,3)][1],
            self.coefficients.sdInfo[(3,3)][0],self.coefficients.sdInfo[(3,3)][1],
            self.coefficients.sdInfo[(3,1)][0],self.coefficients.sdInfo[(3,1)][1],
            self.coefficients.sdInfo[(3,2)][0],self.coefficients.sdInfo[(3,2)][1],
            self.offset[0],self.offset[1],self.offset[2],self.offset[3],
            self.stride[0],self.stride[1],self.stride[2],self.stride[3],
            r,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.coefficients.ebqe_phi,
            self.coefficients.ebqe_n,
            self.coefficients.ebqe_kappa,
            #VRANS start
            self.coefficients.ebqe_porosity,
            #VRANS end
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.isDOFBoundary[1],
            self.numericalFlux.isDOFBoundary[2],
            self.numericalFlux.isDOFBoundary[3],
            self.ebqe[('advectiveFlux_bc_flag',0)],
            self.ebqe[('advectiveFlux_bc_flag',1)],
            self.ebqe[('advectiveFlux_bc_flag',2)],
            self.ebqe[('advectiveFlux_bc_flag',3)],
            self.ebqe[('diffusiveFlux_bc_flag',1,1)],
            self.ebqe[('diffusiveFlux_bc_flag',2,2)],
            self.ebqe[('diffusiveFlux_bc_flag',3,3)],
            self.numericalFlux.ebqe[('u',0)],
            self.ebqe[('advectiveFlux_bc',0)],
            self.ebqe[('advectiveFlux_bc',1)],
            self.ebqe[('advectiveFlux_bc',2)],
            self.ebqe[('advectiveFlux_bc',3)],
            self.numericalFlux.ebqe[('u',1)],
            self.ebqe[('diffusiveFlux_bc',1,1)],
            self.ebqe[('penalty')],
            self.numericalFlux.ebqe[('u',2)],
            self.ebqe[('diffusiveFlux_bc',2,2)],
            self.numericalFlux.ebqe[('u',3)],
            self.ebqe[('diffusiveFlux_bc',3,3)],
            self.q[('velocity',0)],
            self.ebqe[('velocity',0)],
            self.ebq_global[('totalFlux',0)])

	if self.forceStrongConditions:#
	    for cj in range(len(self.dirichletConditionsForceDOF)):#
		for dofN,g in self.dirichletConditionsForceDOF[cj].DOFBoundaryConditionsDict.iteritems():
                     r[self.offset[cj]+self.stride[cj]*dofN] = 0

        if self.stabilization:
            self.stabilization.accumulateSubgridMassHistory(self.q)
        log("Global residual",level=9,data=r)
        self.nonlinear_function_evaluations += 1
    def getJacobian(self,jacobian):
	cfemIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian,
				       jacobian)
        self.rans2p.calculateJacobian(#element
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[1].femSpace.psi,
            self.u[1].femSpace.grad_psi,
            self.u[1].femSpace.psi,
            self.u[1].femSpace.grad_psi,
            #element boundary
            self.u[0].femSpace.elementMaps.psi_trace,
            self.u[0].femSpace.elementMaps.grad_psi_trace,
            self.elementBoundaryQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[1].femSpace.psi_trace,
            self.u[1].femSpace.grad_psi_trace,
            self.u[1].femSpace.psi_trace,
            self.u[1].femSpace.grad_psi_trace,
            self.u[0].femSpace.elementMaps.boundaryNormals,
            self.u[0].femSpace.elementMaps.boundaryJacobians,
            self.mesh.elementDiametersArray,
            self.stabilization.hFactor,
            self.mesh.nElements_global,
            self.coefficients.useRBLES,
            self.timeIntegration.alpha_bdf,
            self.coefficients.epsFact_density,
            self.coefficients.epsFact,
            self.coefficients.sigma,
            self.coefficients.rho_0,
            self.coefficients.nu_0,
            self.coefficients.rho_1,
            self.coefficients.nu_1,
            self.Ct_sge,
            self.Cd_sge,
            self.shockCapturing.shockCapturingFactor,
            #VRANS start
            self.coefficients.linearDragFactor,
            self.coefficients.nonlinearDragFactor,
            self.coefficients.q_porosity,
            self.coefficients.q_meanGrain,
            self.q[('r',0)],
            #VRANS end
            self.u[0].femSpace.dofMap.l2g,
            self.u[1].femSpace.dofMap.l2g,
            self.u[0].dof,
            self.u[1].dof,
            self.u[2].dof,
            self.u[3].dof,
            self.coefficients.g,
            self.coefficients.q_phi,
            self.coefficients.q_n,
            self.coefficients.q_kappa,
            self.timeIntegration.beta_bdf[1],
            self.timeIntegration.beta_bdf[2],
            self.timeIntegration.beta_bdf[3],
            self.stabilization.v_last,
            self.q[('cfl',0)],
            self.shockCapturing.numDiff_last[1],
            self.shockCapturing.numDiff_last[2],
            self.shockCapturing.numDiff_last[3],
            self.coefficients.sdInfo[(1,1)][0],self.coefficients.sdInfo[(1,1)][1],
            self.coefficients.sdInfo[(1,2)][0],self.coefficients.sdInfo[(1,2)][1],
            self.coefficients.sdInfo[(1,3)][0],self.coefficients.sdInfo[(1,3)][1],
            self.coefficients.sdInfo[(2,2)][0],self.coefficients.sdInfo[(2,2)][1],
            self.coefficients.sdInfo[(2,1)][0],self.coefficients.sdInfo[(2,1)][1],
            self.coefficients.sdInfo[(2,3)][0],self.coefficients.sdInfo[(2,3)][1],
            self.coefficients.sdInfo[(3,3)][0],self.coefficients.sdInfo[(3,3)][1],
            self.coefficients.sdInfo[(3,1)][0],self.coefficients.sdInfo[(3,1)][1],
            self.coefficients.sdInfo[(3,2)][0],self.coefficients.sdInfo[(3,2)][1],
            self.csrRowIndeces[(0,0)],self.csrColumnOffsets[(0,0)],
            self.csrRowIndeces[(0,1)],self.csrColumnOffsets[(0,1)],
            self.csrRowIndeces[(0,2)],self.csrColumnOffsets[(0,2)],
            self.csrRowIndeces[(0,3)],self.csrColumnOffsets[(0,3)],
            self.csrRowIndeces[(1,0)],self.csrColumnOffsets[(1,0)],
            self.csrRowIndeces[(1,1)],self.csrColumnOffsets[(1,1)],
            self.csrRowIndeces[(1,2)],self.csrColumnOffsets[(1,2)],
            self.csrRowIndeces[(1,3)],self.csrColumnOffsets[(1,3)],
            self.csrRowIndeces[(2,0)],self.csrColumnOffsets[(2,0)],
            self.csrRowIndeces[(2,1)],self.csrColumnOffsets[(2,1)],
            self.csrRowIndeces[(2,2)],self.csrColumnOffsets[(2,2)],
            self.csrRowIndeces[(2,3)],self.csrColumnOffsets[(2,3)],
            self.csrRowIndeces[(3,0)],self.csrColumnOffsets[(3,0)],
            self.csrRowIndeces[(3,1)],self.csrColumnOffsets[(3,1)],
            self.csrRowIndeces[(3,2)],self.csrColumnOffsets[(3,2)],
            self.csrRowIndeces[(3,3)],self.csrColumnOffsets[(3,3)],
            jacobian,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.coefficients.ebqe_phi,
            self.coefficients.ebqe_n,
            self.coefficients.ebqe_kappa,
            #VRANS start
            self.coefficients.ebqe_porosity,
            #VRANS end
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.isDOFBoundary[1],
            self.numericalFlux.isDOFBoundary[2],
            self.numericalFlux.isDOFBoundary[3],
            self.ebqe[('advectiveFlux_bc_flag',0)],
            self.ebqe[('advectiveFlux_bc_flag',1)],
            self.ebqe[('advectiveFlux_bc_flag',2)],
            self.ebqe[('advectiveFlux_bc_flag',3)],
            self.ebqe[('diffusiveFlux_bc_flag',1,1)],
            self.ebqe[('diffusiveFlux_bc_flag',2,2)],
            self.ebqe[('diffusiveFlux_bc_flag',3,3)],
            self.numericalFlux.ebqe[('u',0)],
            self.ebqe[('advectiveFlux_bc',0)],
            self.ebqe[('advectiveFlux_bc',1)],
            self.ebqe[('advectiveFlux_bc',2)],
            self.ebqe[('advectiveFlux_bc',3)],
            self.numericalFlux.ebqe[('u',1)],
            self.ebqe[('diffusiveFlux_bc',1,1)],
            self.ebqe[('penalty')],
            self.numericalFlux.ebqe[('u',2)],
            self.ebqe[('diffusiveFlux_bc',2,2)],
            self.numericalFlux.ebqe[('u',3)],
            self.ebqe[('diffusiveFlux_bc',3,3)],
            self.csrColumnOffsets_eb[(0,0)],
            self.csrColumnOffsets_eb[(0,1)],
            self.csrColumnOffsets_eb[(0,2)],
            self.csrColumnOffsets_eb[(0,3)],
            self.csrColumnOffsets_eb[(1,0)],
            self.csrColumnOffsets_eb[(1,1)],
            self.csrColumnOffsets_eb[(1,2)],
            self.csrColumnOffsets_eb[(1,3)],
            self.csrColumnOffsets_eb[(2,0)],
            self.csrColumnOffsets_eb[(2,1)],
            self.csrColumnOffsets_eb[(2,2)],
            self.csrColumnOffsets_eb[(2,3)],
            self.csrColumnOffsets_eb[(3,0)],
            self.csrColumnOffsets_eb[(3,1)],
            self.csrColumnOffsets_eb[(3,2)],
            self.csrColumnOffsets_eb[(3,3)])

        #Load the Dirichlet conditions directly into residual
        if self.forceStrongConditions:
            scaling = 1.0#probably want to add some scaling to match non-dirichlet diagonals in linear system 
            for cj in range(self.nc):
                for dofN in self.dirichletConditionsForceDOF[cj].DOFBoundaryConditionsDict.keys():
                    global_dofN = self.offset[cj]+self.stride[cj]*dofN
                    for i in range(self.rowptr[global_dofN],self.rowptr[global_dofN+1]):
                        if (self.colind[i] == global_dofN):
                            #print "RBLES forcing residual cj = %s dofN= %s global_dofN= %s was self.nzval[i]= %s now =%s " % (cj,dofN,global_dofN,self.nzval[i],scaling)
                            self.nzval[i] = scaling
                        else:
                            self.nzval[i] = 0.0
                            #print "RBLES zeroing residual cj = %s dofN= %s global_dofN= %s " % (cj,dofN,global_dofN)
        log("Jacobian ",level=10,data=jacobian)
        #mwf decide if this is reasonable for solver statistics
        self.nonlinear_function_jacobian_evaluations += 1
        return jacobian
    def calculateElementQuadrature(self):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points.
        
        This function should be called only when the mesh changes.
        """
        self.u[0].femSpace.elementMaps.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesRef(self.elementQuadraturePoints)
        self.u[1].femSpace.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[1].femSpace.getBasisGradientValuesRef(self.elementQuadraturePoints)
        self.coefficients.initializeElementQuadrature(self.timeIntegration.t,self.q)
        if self.stabilization != None:
            self.stabilization.initializeElementQuadrature(self.mesh,self.timeIntegration.t,self.q)
            self.stabilization.initializeTimeIntegration(self.timeIntegration)
        if self.shockCapturing != None:
            self.shockCapturing.initializeElementQuadrature(self.mesh,self.timeIntegration.t,self.q)
    def calculateElementBoundaryQuadrature(self):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points on element boundaries.

        This function should be called only when the mesh changes.
        """
        pass
    def calculateExteriorElementBoundaryQuadrature(self):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points on global element boundaries.

        This function should be called only when the mesh changes.
        """
        #
        #get physical locations of element boundary quadrature points
        #
	#assume all components live on the same mesh
        self.u[0].femSpace.elementMaps.getBasisValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[1].femSpace.getBasisValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[1].femSpace.getBasisGradientValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getValuesGlobalExteriorTrace(self.elementBoundaryQuadraturePoints,
                                                                    self.ebqe['x'])
        self.fluxBoundaryConditionsObjectsDict = dict([(cj,FluxBoundaryConditions(self.mesh,
                                                                                  self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                                                  self.ebqe[('x')],
                                                                                  self.advectiveFluxBoundaryConditionsSetterDict[cj],
                                                                                  self.diffusiveFluxBoundaryConditionsSetterDictDict[cj]))
                                                       for cj in self.advectiveFluxBoundaryConditionsSetterDict.keys()])
        self.coefficients.initializeGlobalExteriorElementBoundaryQuadrature(self.timeIntegration.t,self.ebqe)
    def estimate_mt(self):
        pass
    def calculateSolutionAtQuadrature(self):
        pass
    def calculateAuxiliaryQuantitiesAfterStep(self):
        #cek todo memory optimized version
        #self.rans2p.calculateVelocityAverage(self.u[1].femSpace.elementMaps.permutations,
        #     self.mesh.nExteriorElementBoundaries_global,
        #     self.mesh.exteriorElementBoundariesArray,
        #     self.mesh.nInteriorElementBoundaries_global,
        #     self.mesh.interiorElementBoundariesArray,
        #     self.mesh.elementBoundaryElementsArray,
        #     self.mesh.elementBoundaryLocalElementBoundariesArray,
        #     self.u[1].femSpace.dofMap.l2g,
        #     self.u[1].dof,
        #     self.u[2].dof,
        #     self.u[3].dof,
        #     self.ebq[('v',0)], 
        #     self.ebqe[('velocity',0)],
        #     self.ebq_global[('velocityAverage',0)])
        OneLevelTransport.calculateAuxiliaryQuantitiesAfterStep(self)
