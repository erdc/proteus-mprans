from proteus import *
from twp_navier_stokes_p import *
from windtunnel import *

timeIntegration = BackwardEuler_cfl
stepController  = Min_dt_controller

femSpaces = {0:basis,
	     1:basis,
	     2:basis}

massLumping       = False
numericalFluxType = None
conservativeFlux  = None

numericalFluxType = RANS2P.NumericalFlux
subgridError = RANS2P.SubgridError(coefficients,nd,lag=ns_lag_subgridError,hFactor=hFactor)
shockCapturing = RANS2P.ShockCapturing(coefficients,nd,ns_shockCapturingFactor,lag=ns_lag_shockCapturing)

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

tolFac = 0.0
nl_atol_res = 1.0e-4
useEisenstatWalker = True
maxNonlinearIts = 20
maxLineSearches = 0

