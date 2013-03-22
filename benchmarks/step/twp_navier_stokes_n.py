from proteus import *
from twp_navier_stokes_p import *
from step import *

timeIntegration = BackwardEuler_cfl
stepController  = Min_dt_controller

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
        self.penalty_constant = 100.0
        self.includeBoundaryAdjoint=True
        self.boundaryAdjoint_sigma=1.0
        self.hasInterior=False

numericalFluxType = NumericalFluxType
subgridError = NavierStokesASGS_velocity_pressure_optV2(coefficients,nd,lag=True,delayLagSteps=1,hFactor=hFactor,noPressureStabilization=False)
shockCapturing = NavierStokes_SC_opt(coefficients,nd,ns_shockCapturingFactor,lag=True)

fullNewtonFlag = True
multilevelNonlinearSolver = NewtonNS
levelNonlinearSolver      = NewtonNS

nonlinearSmoother = None
linearSmoother    = SimpleNavierStokes3D

matrix = SparseMatrix

if not use_petsc4py:
    multilevelLinearSolver = LU
    levelLinearSolver      = LU
else:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver      = KSP_petsc4py
linear_solver_options_prefix = 'rans2p_'
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest             = 'r-true'

tolFac = 1.0e-3
nl_atol_res = 1.0e-4
useEisenstatWalker = False#True
maxNonlinearIts = 20
maxLineSearches = 0

