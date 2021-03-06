from proteus import *
from proteus.default_n import *
from ls_obstacleInTank_3d_p import *

if useBackwardEuler_ls:
    timeIntegration = BackwardEuler_cfl
    stepController = Min_dt_controller
    stepController = HeuristicNL_dt_controller
    nonlinearIterationsFloor = 3
    nonlinearIterationsCeil=5
    dtNLgrowFactor  = 2.0
    dtNLreduceFactor= 0.5#75
else:
    timeIntegration = FLCBDF
    stepController = FLCBDF_controller_sys
    rtol_u[0] = 1.0e-2
    atol_u[0] = 1.0e-2

if spaceOrder == 1:
    femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
elif spaceOrder == 2:
    femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}

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

tolFac = 1.0e-6#

nl_atol_res = 1.0e-6#0.001*he#1.0e-8#should be linear with lagging

maxNonlinearIts = 50

matrix = SparseMatrix

if usePETSc:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver = KSP_petsc4py
    #multilevelLinearSolver = PETSc
    #levelLinearSolver = PETSc
    linear_solver_options_prefix = 'ncls_'
    #    linearSmoother = StarILU
    linearSmoother = None
    linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU

linTolFac = 0.001

conservativeFlux = None
