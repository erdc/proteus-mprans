from proteus import *
from proteus.default_n import *
from floatingCylinder import *
from vof_floatingCylinder_3d_p import *

elementQuadrature = SimplexGaussQuadrature(nd,quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order)


timeIntegration = BackwardEuler_cfl
stepController=Min_dt_controller
stepController = HeuristicNL_dt_controller
nonlinearIterationsFloor = 2
nonlinearIterationsCeil=2
dtNLgrowFactor  = 2.0
dtNLreduceFactor= 0.5#75

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=vof_shockCapturingFactor,lag=True)#linear
subgridError = Advection_ASGS(coefficients=coefficients,nd=nd,lag=False)
massLumping = False
numericalFluxType = Advection_DiagonalUpwind_IIPG_exterior

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-8

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
