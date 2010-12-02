from proteus import *
from proteus.default_n import *
from obstacleInTank3d import *
from ls_consrv_obstacleInTank_3d_p import *


stepController = FixedStep
timeIntegration = NoIntegration

if useHex:
	femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}

	elementQuadrature = CubeGaussQuadrature(nd,obstacleInTank_quad_order)
	elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,obstacleInTank_quad_order)
else:
	femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

	elementQuadrature = SimplexGaussQuadrature(nd,obstacleInTank_quad_order)
	elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,obstacleInTank_quad_order)

subgridError = None

massLumping = False

numericalFluxType = DoNothing

shockCapturing = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.001

nl_atol_res = 0.0

maxNonlinearIts = 10
maxLineSearches =0

matrix = SparseMatrix

if usePETSc:
    multilevelLinearSolver = PETSc #KSP_petsc4py
    levelLinearSolver = PETSc #KSP_petsc4py
    linear_solver_options_prefix = 'mcorr_'
    linearSmoother = None
    #linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU


linTolFac = 1.0e-6

conservativeFlux = None
