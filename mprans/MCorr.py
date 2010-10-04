import proteus
from proteus.mprans.cMCorr import *

class Coefficients(proteus.TransportCoefficients.TC_base):
    from proteus.ctransportCoefficients import levelSetConservationCoefficientsEvaluate
    from proteus.ctransportCoefficients import levelSetConservationCoefficientsEvaluate_sd
    def __init__(self,applyCorrection=True,epsFactHeaviside=0.0,epsFactDirac=1.0,epsFactDiffusion=2.0,LSModel_index=3,V_model=2,me_model=5,VOFModel_index=4,checkMass=True,sd=True,nd=None,applyCorrectionToDOF=True):
        self.sd=sd
        self.checkMass=checkMass
        self.variableNames=['phiCorr']
        nc=1
        mass={}
        advection={}
        hamiltonian={}
        diffusion={0:{0:{0:'constant'}}}
        potential={0:{0:'u'}}
        reaction={0:{0:'nonlinear'}}
        #reaction={}
        if self.sd:
            assert nd!=None,"You must set the number of dimensions to use sparse diffusion in LevelSetConservationCoefficients"
            sdInfo = {(0,0):(numpy.arange(start=0,stop=nd+1,step=1,dtype='i'),
                             numpy.arange(start=0,stop=nd,step=1,dtype='i'))}
        else:
            sdInfo={}
        TC_base.__init__(self,
                         nc,
                         mass,
                         advection,
                         diffusion,
                         potential,
                         reaction,
                         hamiltonian,
                         self.variableNames,
                         sparseDiffusionTensors=sdInfo,
                         useSparseDiffusion = sd)
        self.levelSetModelIndex=LSModel_index
        self.flowModelIndex=V_model
        self.epsFactHeaviside=epsFactHeaviside
        self.epsFactDirac=epsFactDirac
        self.epsFactDiffusion=epsFactDiffusion
        self.me_model=me_model
        self.VOFModelIndex=VOFModel_index
        self.useC = True
        self.applyCorrection=applyCorrection
        if self.applyCorrection:
            self.applyCorrectionToDOF=applyCorrectionToDOF
        else:
            self.applyCorrection = False
    def initializeMesh(self,mesh):
        self.h=mesh.h
        self.epsHeaviside = self.epsFactHeaviside*mesh.h
        self.epsDirac = self.epsFactDirac*mesh.h
        self.epsDiffusion = self.epsFactDiffusion*mesh.h
    def attachModels(self,modelList):
        import copy
        log("Attaching models in LevelSetConservation")
        #level set
        self.lsModel = modelList[self.levelSetModelIndex]
        self.q_u_ls    = modelList[self.levelSetModelIndex].q[('u',0)]
        self.ebqe_u_ls = modelList[self.levelSetModelIndex].ebqe[('u',0)]
        if modelList[self.levelSetModelIndex].ebq.has_key(('u',0)):
            self.ebq_u_ls = modelList[self.levelSetModelIndex].ebq[('u',0)]
        else:
            self.ebq_u_ls = None
        #volume of fluid
        self.vofModel = modelList[self.VOFModelIndex]
        self.q_H_vof = modelList[self.VOFModelIndex].q[('u',0)]
        self.ebqe_H_vof = modelList[self.VOFModelIndex].ebqe[('u',0)]
        if modelList[self.VOFModelIndex].ebq.has_key(('u',0)):
            self.ebq_H_vof = modelList[self.VOFModelIndex].ebq[('u',0)]
        else:
            self.ebq_H_vof = None
        #correction
        self.massCorrModel = modelList[self.me_model]
        if self.checkMass:
            self.m_tmp = copy.deepcopy(self.massCorrModel.q[('r',0)])
            if self.checkMass:
                self.vofGlobalMass = Norms.scalarDomainIntegral(self.vofModel.q['dV'],
                                                                self.vofModel.q[('u',0)],
                                                                self.massCorrModel.mesh.nElements_owned)
                self.lsGlobalMass = Norms.scalarHeavisideDomainIntegral(self.vofModel.q['dV'],
                                                                        self.lsModel.q[('u',0)],
                                                                        self.massCorrModel.mesh.nElements_owned)
                log("Attach Models MCorr: mass correction %12.5e" % (Norms.scalarDomainIntegral(self.vofModel.q['dV'],
                                                                                                self.massCorrModel.q[('r',0)],
                                                                                                self.massCorrModel.mesh.nElements_owned),),level=2)
                self.fluxGlobal = 0.0
                self.totalFluxGlobal = 0.0
                self.vofGlobalMassArray = [self.vofGlobalMass]
                self.lsGlobalMassArray = [self.lsGlobalMass]
                self.vofGlobalMassErrorArray = [self.vofGlobalMass - self.vofGlobalMassArray[0] + self.vofModel.timeIntegration.dt*self.vofModel.coefficients.fluxIntegral]
                self.lsGlobalMassErrorArray = [self.lsGlobalMass - self.lsGlobalMassArray[0] + self.vofModel.timeIntegration.dt*self.vofModel.coefficients.fluxIntegral]
                self.fluxArray = [self.vofModel.coefficients.fluxIntegral]
                self.timeArray = [self.vofModel.timeIntegration.t]
                log("Attach Models MCorr: Phase 0 mass after mass correction (VOF) %12.5e" % (self.vofGlobalMass,),level=2)
                log("Attach Models MCorr: Phase 0 mass after mass correction (LS) %12.5e" % (self.lsGlobalMass,),level=2)
                log("Attach Models MCorr: Phase  0 mass conservation (VOF) after step = %12.5e" % (self.vofGlobalMass - self.vofModel.coefficients.m_pre + self.vofModel.timeIntegration.dt*self.vofModel.coefficients.fluxIntegral,),level=2)
                log("Attach Models MCorr: Phase  0 mass conservation (LS) after step = %12.5e" % (self.lsGlobalMass - self.lsModel.coefficients.m_pre + self.vofModel.timeIntegration.dt*self.vofModel.coefficients.fluxIntegral,),level=2)
    def initializeElementQuadrature(self,t,cq):
        if self.sd and cq.has_key(('a',0,0)):
            cq[('a',0,0)].fill(self.epsDiffusion)
    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        if self.sd and cebq.has_key(('a',0,0)):
            cebq[('a',0,0)].fill(self.epsDiffusion)
    def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
        if self.sd and cebqe.has_key(('a',0,0)):
            cebqe[('a',0,0)].fill(self.epsDiffusion)
    def preStep(self,t,firstStep=False):
        if self.checkMass:
            log("Phase 0 mass before mass correction (VOF) %12.5e" % (Norms.scalarDomainIntegral(self.vofModel.q['dV'],
                                                                                                 self.vofModel.q[('m',0)],
                                                                                                 self.massCorrModel.mesh.nElements_owned),),level=2)
            log("Phase 0 mass before mass correction (LS) %12.5e" % (Norms.scalarHeavisideDomainIntegral(self.vofModel.q['dV'],
                                                                                                         self.lsModel.q[('m',0)],
                                                                                                         self.massCorrModel.mesh.nElements_owned),),level=2)
        copyInstructions = {'clear_uList':True}
        return copyInstructions
    def postStep(self,t,firstStep=False):
        if self.applyCorrection:
            
            self.vofModel.q[('m',0)] += self.massCorrModel.q[('r',0)]
            self.lsModel.q[('m',0)] += self.massCorrModel.q[('u',0)]
            if self.vofModel.q[('u',0)] is not self.vofModel.q[('m',0)]:
                self.vofModel.q[('u',0)][:]=self.vofModel.q[('m',0)]
            if self.lsModel.q[('u',0)] is not self.lsModel.q[('m',0)]:
                self.lsModel.q[('u',0)][:]=self.lsModel.q[('m',0)]
            if self.vofModel.q.has_key(('mt',0)):
                self.vofModel.timeIntegration.calculateElementCoefficients(self.vofModel.q)
                self.vofModel.timeIntegration.lastStepErrorOk()
            if self.applyCorrectionToDOF:
                self.lsModel.u[0].dof += self.massCorrModel.u[0].dof
            if self.lsModel.q.has_key(('mt',0)):
                self.lsModel.timeIntegration.calculateElementCoefficients(self.lsModel.q)
                self.lsModel.timeIntegration.lastStepErrorOk()
            if self.checkMass:
                self.vofGlobalMass = Norms.scalarDomainIntegral(self.vofModel.q['dV'],
                                                                self.vofModel.q[('u',0)],
                                                                self.massCorrModel.mesh.nElements_owned)
                self.lsGlobalMass = Norms.scalarHeavisideDomainIntegral(self.vofModel.q['dV'],
                                                                        self.lsModel.q[('u',0)],
                                                                        self.massCorrModel.mesh.nElements_owned)
                log("mass correction %12.5e" % (Norms.scalarDomainIntegral(self.vofModel.q['dV'],
                                                                           self.massCorrModel.q[('r',0)],
                                                                           self.massCorrModel.mesh.nElements_owned),),level=2)
                self.fluxGlobal = self.vofModel.coefficients.fluxIntegral*self.vofModel.timeIntegration.dt
                self.totalFluxGlobal += self.vofModel.coefficients.fluxIntegral*self.vofModel.timeIntegration.dt
                self.vofGlobalMassArray.append(self.vofGlobalMass)
                self.lsGlobalMassArray.append(self.lsGlobalMass)
                self.vofGlobalMassErrorArray.append(self.vofGlobalMass - self.vofGlobalMassArray[0] + self.totalFluxGlobal)
                self.lsGlobalMassErrorArray.append(self.lsGlobalMass - self.lsGlobalMassArray[0] + self.totalFluxGlobal)
                self.fluxArray.append(self.vofModel.coefficients.fluxIntegral)
                self.timeArray.append(self.vofModel.timeIntegration.t)
                log("Phase 0 mass after mass correction (VOF) %12.5e" % (self.vofGlobalMass,),level=2)
                log("Phase 0 mass after mass correction (LS) %12.5e" % (self.lsGlobalMass,),level=2)
                log("Phase  0 mass conservation (VOF) after step = %12.5e" % (self.vofGlobalMass - self.vofModel.coefficients.m_last + self.fluxGlobal,),level=2)
                log("Phase  0 mass conservation (LS) after step = %12.5e" % (self.lsGlobalMass - self.lsModel.coefficients.m_last + self.fluxGlobal,),level=2)
        copyInstructions = {}
        return copyInstructions
    def evaluate(self,t,c):
        import math
        if c[('u',0)].shape == self.q_u_ls.shape:
            u_ls = self.q_u_ls
            H_vof = self.q_H_vof
        elif c[('u',0)].shape == self.ebqe_u_ls.shape:
            u_ls = self.ebqe_u_ls
            H_vof = self.ebqe_H_vof
        elif self.ebq_u_ls != None and c[('u',0)].shape == self.ebq_u_ls.shape:
            u_ls = self.ebq_u_ls
            H_vof = self.ebq_H_vof
        else:
            #\todo trap errors in TransportCoefficients.py
            u_ls = None
            H_vof = None
        if u_ls != None and H_vof != None:
            if self.useC:
                if self.sd:
                    self.levelSetConservationCoefficientsEvaluate_sd(self.epsHeaviside,
                                                                     self.epsDirac,
                                                                     u_ls,
                                                                     H_vof,
                                                                     c[('u',0)],
                                                                     c[('r',0)],
                                                                     c[('dr',0,0)])
                else:
                    self.levelSetConservationCoefficientsEvaluate(self.epsHeaviside,
                                                                  self.epsDirac,
                                                                  self.epsDiffusion,
                                                                  u_ls,
                                                                  H_vof,
                                                                  c[('u',0)],
                                                                  c[('r',0)],
                                                                  c[('dr',0,0)],
                                                                  c[('a',0,0)])
        if (self.checkMass and c[('u',0)].shape == self.q_u_ls.shape):
            self.m_tmp[:] = H_vof
            self.m_tmp += self.massCorrModel.q[('r',0)]
            log("mass correction during Newton %12.5e" % (Norms.scalarDomainIntegral(self.vofModel.q['dV'],
                                                                                     self.massCorrModel.q[('r',0)],
                                                                                     self.massCorrModel.mesh.nElements_owned),),level=2)
            log("Phase 0 mass during Newton %12.5e" % (Norms.scalarDomainIntegral(self.vofModel.q['dV'],
                                                                                 self.m_tmp,
                                                                                  self.massCorrModel.mesh.nElements_owned),),level=2)

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
                 stressTraceBoundaryConditionsSetterDict=None,
                 stabilization=None,
                 shockCapturing=None,
                 conservativeFluxDict=None,
                 numericalFluxType=None,
                 TimeIntegrationClass=None,
                 massLumping=False,
                 reactionLumping=False,
                 options=None,
                 name='defaultName',
                 reuse_trial_and_test_quadrature=True,
                 sd = True,
                 movingDomain=False):#,
        from proteus import Comm
        #
        #set the objects describing the method and boundary conditions
        #
        self.movingDomain=movingDomain
        self.tLast_mesh=None
        #
        self.name=name
        self.sd=sd
        self.Hess=False
        self.lowmem=True
        self.timeTerm=True#allow turning off  the  time derivative
        #self.lowmem=False
        self.testIsTrial=True
        self.phiTrialIsTrial=True            
        self.u = uDict
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
        if type(self.u[0].femSpace) == C0_AffineLinearOnSimplexWithNodalBasis:
            print self.nQuadraturePoints_element
            if self.nSpace_global == 3:
                assert(self.nQuadraturePoints_element == 5)
            elif self.nSpace_global == 2:
                assert(self.nQuadraturePoints_element == 6)
            elif self.nSpace_global == 1:
                assert(self.nQuadraturePoints_element == 3)

            print self.nElementBoundaryQuadraturePoints_elementBoundary
            if self.nSpace_global == 3:
                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 4)
            elif self.nSpace_global == 2:
                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 4)
            elif self.nSpace_global == 1:
                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 1)

        #pdb.set_trace()
        #
        #simplified allocations for test==trial and also check if space is mixed or not
        #
        self.q={}
        self.ebq={}
        self.ebq_global={}
        self.ebqe={}
        self.phi_ip={}
        #mesh
        self.q[('u',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('grad(u)',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),'d')
        self.q[('r',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.points_elementBoundaryQuadrature= set()
        self.scalars_elementBoundaryQuadrature= set([('u',ci) for ci in range(self.nc)])
        self.vectors_elementBoundaryQuadrature= set()
        self.tensors_elementBoundaryQuadrature= set()
        log(memory("element and element boundary Jacobians","OneLevelTransport"),level=4)
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
        self.globalResidualDummy = None
        compKernelFlag = 0
        self.mcorr = cMCorr_base(self.nSpace_global,
                                 self.nQuadraturePoints_element,
                                 self.u[0].femSpace.elementMaps.localFunctionSpace.dim,
                                 self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.dim,
                                 self.testSpace[0].referenceFiniteElement.localFunctionSpace.dim,
                                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                                 compKernelFlag)
    #mwf these are getting called by redistancing classes,
    def calculateCoefficients(self):
        pass
    def calculateElementResidual(self):
        if self.globalResidualDummy != None:
            self.getResidual(self.u[0].dof,self.globalResidualDummy)
    def getResidual(self,u,r):
        import pdb
        import copy
        """
        Calculate the element residuals and add in to the global residual
        """
        r.fill(0.0)
        #Load the unknowns into the finite element dof
        self.setUnknowns(u)
        #no flux boundary conditions
        self.mcorr.calculateResidual(#element
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            #element boundary
            self.u[0].femSpace.elementMaps.psi_trace,
            self.u[0].femSpace.elementMaps.grad_psi_trace,
            self.elementBoundaryQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.elementMaps.boundaryNormals,
            self.u[0].femSpace.elementMaps.boundaryJacobians,
            #physics
            self.mesh.nElements_global,
                                 self.coefficients.epsFactHeaviside,
                                 self.coefficients.epsFactDirac,
                                 self.coefficients.epsFactDiffusion,
                                 self.u[0].femSpace.dofMap.l2g,
                                 self.mesh.elementDiametersArray,
                                 self.u[0].dof,
                                 self.coefficients.q_u_ls,
                                 self.coefficients.q_H_vof,
                                 self.q[('u',0)],
                                 self.q[('r',0)],
                                 self.offset[0],self.stride[0],
                                 r)
        log("Global residual",level=9,data=r)
        self.nonlinear_function_evaluations += 1
        if self.globalResidualDummy == None:
            self.globalResidualDummy = numpy.zeros(r.shape,'d')
    def getJacobian(self,jacobian):
        #import superluWrappers
        #import numpy
        import pdb
	cfemIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian,jacobian)
        self.mcorr.calculateJacobian(#element
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            #element boundary
            self.u[0].femSpace.elementMaps.psi_trace,
            self.u[0].femSpace.elementMaps.grad_psi_trace,
            self.elementBoundaryQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.elementMaps.boundaryNormals,
            self.u[0].femSpace.elementMaps.boundaryJacobians,
            self.mesh.nElements_global,
                                 self.coefficients.epsFactHeaviside,
                                 self.coefficients.epsFactDirac,
                                 self.coefficients.epsFactDiffusion,
                                 self.u[0].femSpace.dofMap.l2g,
                                 self.mesh.elementDiametersArray,
                                 self.u[0].dof,
                                 self.coefficients.q_u_ls,
                                 self.coefficients.q_H_vof,
                                 self.csrRowIndeces[(0,0)],self.csrColumnOffsets[(0,0)],
                                 jacobian)
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
        self.coefficients.initializeElementQuadrature(self.timeIntegration.t,self.q)
        if self.stabilization != None:
            self.stabilization.initializeElementQuadrature(self.mesh,self.timeIntegration.t,self.q)
            self.stabilization.initializeTimeIntegration(self.timeIntegration)
        if self.shockCapturing != None:
            self.shockCapturing.initializeElementQuadrature(self.mesh,self.timeIntegration.t,self.q)
    def calculateElementBoundaryQuadrature(self):
        pass
    def calculateExteriorElementBoundaryQuadrature(self):
        self.u[0].femSpace.elementMaps.getBasisValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesTraceRef(self.elementBoundaryQuadraturePoints)
    def estimate_mt(self):
        pass
    def calculateSolutionAtQuadrature(self):
        pass
    def calculateAuxiliaryQuantitiesAfterStep(self):
        pass
