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
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest         = 'r-true'

tolFac = 0.0
linTolFac = 0.0
nl_atol_res = he**2
l_atol_res = 0.001*nl_atol_res
useEisenstatWalker = True

maxNonlinearIts = 50
maxLineSearches = 0

