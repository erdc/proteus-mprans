from proteus import *
from ls_p import *

timeIntegration = BackwardEuler
stepController  = Min_dt_controller

femSpaces = {0:basis}

massLumping       = False
numericalFluxType = None
conservativeFlux  = None
numericalFluxType = NCLS.NumericalFlux
subgridError      = NCLS.SubgridError(coefficients,nd)
shockCapturing    = NCLS.ShockCapturing(coefficients,nd,shockCapturingFactor=ls_shockCapturingFactor,lag=ls_lag_shockCapturing)

fullNewtonFlag  = True
multilevelNonlinearSolver = Newton
levelNonlinearSolver      = Newton

nonlinearSmoother = None
linearSmoother    = None

matrix = SparseMatrix
if not use_petsc4py:
    multilevelLinearSolver = LU
    levelLinearSolver      = LU
else:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver      = KSP_petsc4py
linear_solver_options_prefix = 'ncls_'
levelNonlinearSolverConvergenceTest = 'rits'
linearSolverConvergenceTest         = 'rits'

tolFac = 1e-3
nl_atol_res = 0.0

maxNonlinearIts = 1
maxLineSearches = 0

