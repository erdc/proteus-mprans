from proteus import *
from proteus.default_n import *
from sloshbox3d import *
from vof_sloshbox_3d_p import *



timeIntegration = BackwardEuler_cfl
stepController=Min_dt_controller
if timeOrder == 2:
    timeIntegration = VBDF
    stepController = Min_dt_cfl_controller
#timeIntegration = BackwardEuler
#stepController=FixedStep
if useHex:
    if spaceOrder == 1:
        femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}
    elif spaceOrder == 2:
        femSpaces = {0:C0_AffineLagrangeOnCubeWithNodalBasis}
    elementQuadrature = CubeGaussQuadrature(nd,sloshbox_quad_order)
    elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,sloshbox_quad_order)
else:
    if spaceOrder == 1:
        femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
    elif spaceOrder == 2:
        femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
    elementQuadrature = SimplexGaussQuadrature(nd,sloshbox_quad_order)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,sloshbox_quad_order)

shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=vof_shockCapturingFactor,lag=False)#linear
subgridError = Advection_ASGS(coefficients=coefficients,nd=nd,lag=False)
massLumping = False
numericalFluxType = Advection_DiagonalUpwind_IIPG_exterior
conservativeFlux = None


multilevelNonlinearSolver  = Newton
levelNonlinearSolver = Newton

fullNewtonFlag = True

nonlinearSmoother = None
linearSmoother = None

tolFac      = 1e-4
nl_atol_res = 0.0

levelNonlinearSolverConvergenceTest='rits'
maxNonlinearIts = 10
maxLineSearches = 0

matrix = SparseMatrix

multilevelLinearSolver = PETSc
levelLinearSolver = PETSc
linear_solver_options_prefix = 'vof_'

levelNonlinearSolverConvergenceTest='rits'
linearSolverConvergenceTest = 'rits'

