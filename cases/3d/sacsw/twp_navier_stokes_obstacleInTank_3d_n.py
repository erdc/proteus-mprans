from proteus import *
from proteus.default_n import *
from twp_navier_stokes_obstacleInTank_3d_p import *
from obstacleInTank3d import *

if useBackwardEuler:
    timeIntegration = BackwardEuler_cfl
    stepController = Min_dt_controller
    stepController = HeuristicNL_dt_controller
    nonlinearIterationsFloor = 2
    nonlinearIterationsCeil=4
    nonlinearIterationsFloor = 3
    nonlinearIterationsCeil=4
    dtNLgrowFactor  = 2.0
    dtNLreduceFactor= 0.5#75
else:
    timeIntegration = FLCBDF
    stepController = FLCBDF_controller_sys
    rtol_u[1] = 1.0e-2
    rtol_u[2] = 1.0e-2
    rtol_u[3] = 1.0e-2
    atol_u[1] = 1.0e-2
    atol_u[2] = 1.0e-2
    atol_u[3] = 1.0e-2

noPressureStabilization=False
if spaceOrder==1:
    femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
                 1:C0_AffineLinearOnSimplexWithNodalBasis,
                 2:C0_AffineLinearOnSimplexWithNodalBasis,
                 3:C0_AffineLinearOnSimplexWithNodalBasis}
    hFactor=1.0
if spaceOrder==2:
    femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis,
                 1:C0_AffineQuadraticOnSimplexWithNodalBasis,
                 2:C0_AffineQuadraticOnSimplexWithNodalBasis,
                 3:C0_AffineQuadraticOnSimplexWithNodalBasis}
    hFactor=0.5
elementQuadrature = SimplexGaussQuadrature(nd,obstacleInTank_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,obstacleInTank_quad_order)

subgridError = None

subgridError = NavierStokesASGS_velocity_pressure_optV2(coefficients,nd,lag=lag_ns_subgridError,delayLagSteps=1,hFactor=hFactor,noPressureStabilization=noPressureStabilization)

massLumping = False

shockCapturing = NavierStokes_SC_opt(coefficients,nd,ns_shockCapturingFactor,lag=lag_ns_shockCapturing)

numericalFluxType = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

maxNonlinearIts = 10
maxLineSearches =0

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-4#0.0001*he

matrix = SparseMatrix

numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior #need weak for parallel and global conservation

if usePETSc:    
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver = KSP_petsc4py
    linear_solver_options_prefix = 'rans2p_'
#    linearSmoother = StarILU
    linearSmoother = None
    linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU

linTolFac = 0.001

#conservativeFlux = {0:'pwl-bdm-opt'}
#conservativeFlux = {0:'pwl'}#,1:'pwl-bdm',2:'pwl-bdm'}
#conservativeFlux = {0:'pwc'}#,1:'pwc',2:'pwc'}
#conservativeFlux = {0:'point-eval',1:'point-eval',2:'point-eval'}
#auxiliaryVariables=[rc]
