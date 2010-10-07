from proteus import *
from proteus.default_n import *
from obstacleInTank3d import *
from vof_obstacleInTank_3d_p import *

elementQuadrature = SimplexGaussQuadrature(nd,obstacleInTank_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,obstacleInTank_quad_order)


timeIntegration = BackwardEuler_cfl
stepController=Min_dt_controller
stepController = HeuristicNL_dt_controller
nonlinearIterationsFloor = 3
nonlinearIterationsCeil=5
dtNLgrowFactor  = 2.0
dtNLreduceFactor= 0.5#75
#stepController = HeuristicNL_dt_controller
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=vof_shockCapturingFactor,lag=lag_vof_shockCapturing)#linear
subgridError = Advection_ASGS(coefficients=coefficients,nd=nd,lag=False)#it's linear anyway
massLumping = False
numericalFluxType = Advection_DiagonalUpwind_IIPG_exterior

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 0.001*he#1.0e-8

maxNonlinearIts = 50

matrix = SparseMatrix

if usePETSc:
     multilevelLinearSolver = KSP_petsc4py
     levelLinearSolver = KSP_petsc4py
     linear_solver_options_prefix = 'vof_'
#     linearSmoother = StarILU
     linearSmoother = None
     linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU

linTolFac = 0.001

conservativeFlux = None
