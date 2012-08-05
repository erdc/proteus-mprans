from proteus import *
from wavetank import *
from ls_consrv_p import *

timeIntegrator  = ForwardIntegrator
timeIntegration = NoIntegration

femSpaces = {0:basis}

subgridError      = None
massLumping       = False
numericalFluxType = DoNothing
conservativeFlux  = None
shockCapturing    = None

fullNewtonFlag = True
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

linear_solver_options_prefix = 'mcorr_'
linearSolverConvergenceTest  = 'rits'

tolFac = 0.0
nl_atol_res = 1.0e-5

maxNonlinearIts = 10
maxLineSearches = 0
