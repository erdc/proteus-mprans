from proteus import *
from kappa_p import *

timeIntegration = BackwardEuler_cfl
stepController  = Min_dt_cfl_controller

femSpaces = {0:basis}

massLumping       = False
numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior
conservativeFlux  = None
subgridError      = Advection_ASGS(coefficients,nd,lag=False) #needs to be addressed or just skip because going to handle in optimized code anyway?
shockCapturing    = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=kappa_shockCapturingFactor,lag=True)

fullNewtonFlag  = True
multilevelNonlinearSolver = Newton
levelNonlinearSolver      = Newton

nonlinearSmoother = None
linearSmoother    = None
#printNonlinearSolverInfo = True
matrix = SparseMatrix
if use_petsc4py:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver      = KSP_petsc4py
else:
    multilevelLinearSolver = LU
    levelLinearSolver      = LU

linear_solver_options_prefix = 'kappa_'
levelNonlinearSolverConvergenceTest = 'r'#'rits'
linearSolverConvergenceTest         = 'r'#'rits'

tolFac = 0.0
nl_atol_res = 1.0e-6
nl_rtol_res = 0.0

maxNonlinearIts = 10
maxLineSearches = 0

