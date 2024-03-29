from proteus import *
from proteus.default_n import *
from obstacleInTank3d import *
from vof_obstacleInTank_3d_p import *

elementQuadrature = CubeGaussQuadrature(nd,obstacleInTank_quad_order)

elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,obstacleInTank_quad_order)


timeIntegration = BackwardEuler_cfl
stepController = FixedStep

if spaceOrder == 1:
    femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}
elif spaceOrder == 2:
    femSpaces = {0:C0_AffineLagrangeOnCubeWithNodalBasis}


shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=vof_shockCapturingFactor,lag=lag_vof_shockCapturing)#linear
subgridError = Advection_ASGS(coefficients=coefficients,nd=nd,lag=False)#it's linear anyway
massLumping = False
numericalFluxType = Advection_DiagonalUpwind_IIPG_exterior

multilevelNonlinearSolver  = Newton #NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0001

nl_atol_res = 0.0

maxNonlinearIts = 10

matrix = SparseMatrix

nonlinearSolverConvergenceTest = 'rits'
levelNonlinearSolverConvergenceTest = 'rits'

if usePETSc:
     multilevelLinearSolver = PETSc #KSP_petsc4py
     levelLinearSolver = PETSc #KSP_petsc4py
     linear_solver_options_prefix = 'vof_'
     #linearSmoother = None
     linearSmoother = None ##StarILU
     #linearSolverConvergenceTest = 'r-true'
     nonlinearSolverConvergenceTest='rits'
     levelNonlinearSolverConvergenceTest='rits'
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU

linTolFac = 0.001

conservativeFlux = None
