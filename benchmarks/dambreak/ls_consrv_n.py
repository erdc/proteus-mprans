from proteus import *
from dambreak import *
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

multilevelLinearSolver = PETSc
levelLinearSolver      = PETSc
linear_solver_options_prefix = 'mcorr_'
linearSolverConvergenceTest  = 'rits'

tolFac = 0.0
nl_atol_res = 1.0e-5

maxNonlinearIts = 10
maxLineSearches = 0
