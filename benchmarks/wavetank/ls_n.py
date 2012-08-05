from proteus import *
from ls_p import *

timeIntegration = BackwardEuler
stepController  = Min_dt_controller

femSpaces = {0:basis}

massLumping       = False
numericalFluxType = None
conservativeFlux  = None
numericalFluxType = DoNothing
subgridError      = HamiltonJacobi_ASGS_opt(coefficients,nd,lag=False)
shockCapturing    = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=ls_shockCapturingFactor,lag=True)

fullNewtonFlag  = True
multilevelNonlinearSolver = Newton
levelNonlinearSolver      = Newton

nonlinearSmoother = None
linearSmoother    = None

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

linear_solver_options_prefix = 'ncls_'
levelNonlinearSolverConvergenceTest = 'rits'
linearSolverConvergenceTest         = 'rits'

tolFac = 1e-3
nl_atol_res = 0.0

maxNonlinearIts = 10
maxLineSearches = 0

