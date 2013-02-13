from proteus import *
from twp_navier_stokes_p import *
from dtmb import *

timeIntegration = BackwardEuler
stepController  = Min_dt_controller

timeIntegration = BackwardEuler_cfl
stepController = Min_dt_controller

femSpaces = {0:basis,
	     1:basis,
	     2:basis,
	     3:basis}

massLumping       = False
numericalFluxType = None
conservativeFlux  = None

class NumericalFluxType(NavierStokes_Advection_DiagonalUpwind_Diffusion_SIPG_exterior):
    hasInterior=False
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None):
        NavierStokes_Advection_DiagonalUpwind_Diffusion_SIPG_exterior.__init__(self,vt,getPointwiseBoundaryConditions,
                                                                               getAdvectiveFluxBoundaryConditions,
                                                                               getDiffusiveFluxBoundaryConditions,getPeriodicBoundaryConditions)
        self.penalty_constant = 2.0
        self.includeBoundaryAdjoint=True
        self.boundaryAdjoint_sigma=1.0
        self.hasInterior=False

numericalFluxType = NumericalFluxType
subgridError = NavierStokesASGS_velocity_pressure_optV2(coefficients,nd,lag=True,delayLagSteps=1,hFactor=hFactor,noPressureStabilization=False)
shockCapturing = NavierStokes_SC_opt(coefficients,nd,ns_shockCapturingFactor,lag=True)

fullNewtonFlag = True
multilevelNonlinearSolver = NewtonNS
levelNonlinearSolver      = NewtonNS
multilevelNonlinearSolver = Newton
levelNonlinearSolver      = Newton

nonlinearSmoother = None
linearSmoother    = SimpleNavierStokes3D

matrix = SparseMatrix

if useOldPETSc:
    multilevelLinearSolver = PETSc
    levelLinearSolver      = PETSc
else:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver      = KSP_petsc4py

if useSuperlu:
    multilevelLinearSolver = LU
    levelLinearSolver      = LU
    
linear_solver_options_prefix = 'rans2p_'
levelNonlinearSolverConvergenceTest = 'rits'
linearSolverConvergenceTest         = 'r-true'

tolFac = 0.0
nl_atol_res = 1.0e-3

maxNonlinearIts = 20
maxLineSearches = 0

#auxiliaryVariables = [RelaxationZoneWaveGenerator(twpflowVelocity_w,twpflowVelocity_w,twpflowVelocity_w,xRelaxCenter)]
