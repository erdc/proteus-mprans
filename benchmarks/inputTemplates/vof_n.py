from proteus import *
from example import *
from vof_p import *

timeIntegration = BackwardEuler
timeIntegration = BackwardEuler_cfl
stepController  = Min_dt_controller

femSpaces = {0:basis}

massLumping       = False
numericalFluxType = Advection_DiagonalUpwind_IIPG_exterior
conservativeFlux  = None
shockCapturing    = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=vof_shockCapturingFactor,lag=True)
subgridError      = Advection_ASGS(coefficients=coefficients,nd=nd,lag=False)

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

linear_solver_options_prefix = 'vof_'
levelNonlinearSolverConvergenceTest = 'rits'
linearSolverConvergenceTest         = 'r-true'

tolFac      = 0.0
nl_atol_res = 1.0e-3

maxNonlinearIts = 10
maxLineSearches = 0
