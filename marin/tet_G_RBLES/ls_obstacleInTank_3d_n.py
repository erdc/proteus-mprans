from proteus import *
from proteus.default_n import *
from ls_obstacleInTank_3d_p import *

if useBackwardEuler_ls:
    timeIntegration = BackwardEuler_cfl
    stepController = FixedStep
else:
    timeIntegration = FLCBDF
    stepController = FLCBDF_controller_sys
    rtol_u[0] = 1.0e-2
    atol_u[0] = 1.0e-2

if useHex:
	femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}

	elementQuadrature = CubeGaussQuadrature(nd,obstacleInTank_quad_order)
	elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,obstacleInTank_quad_order)
else:
	femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

	elementQuadrature = SimplexGaussQuadrature(nd,obstacleInTank_quad_order)
	elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,obstacleInTank_quad_order)


subgridError = HamiltonJacobi_ASGS_opt(coefficients,nd,lag=False)#it's  linear anyway

massLumping = False

numericalFluxType = DoNothing

shockCapturing = None

shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=ls_shockCapturingFactor,lag=lag_ls_shockCapturing)

multilevelNonlinearSolver  = Newton#NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0001

nl_atol_res = 0.0#1.0e-8#should be linear with lagging

maxNonlinearIts = 50

matrix = SparseMatrix

nonlinearSolverConvergenceTest = 'rits'
levelNonlinearSolverConvergenceTest = 'rits'

if usePETSc:
    multilevelLinearSolver = PETSc #KSP_petsc4py
    levelLinearSolver = PETSc #KSP_petsc4py
    linear_solver_options_prefix = 'ncls_'
    linearSmoother = None
    #linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU

linTolFac = 0.001

conservativeFlux = None
