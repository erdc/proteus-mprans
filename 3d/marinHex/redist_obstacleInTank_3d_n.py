from proteus import *
from proteus.default_n import *
from redist_obstacleInTank_3d_p import *
from obstacleInTank3d import *

if rdtimeIntegration == 'newton':    
    timeIntegration = NoIntegration
    stepController = Newton_controller
else:
    timeIntegration = BackwardEuler_cfl
    stepController = Osher_PsiTC_controller2
    runCFL=1.0
    rtol_res[0] = 0.001
    atol_res[0] = 0.0

if spaceOrder == 1:
    femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}
if spaceOrder == 2:
    femSpaces = {0:C0_AffineQuadraticOnCubeWithNodalBasis}

elementQuadrature = CubeGaussQuadrature(nd,obstacleInTank_quad_order)

elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,obstacleInTank_quad_order)

if rdtimeIntegration != 'newton':    
    subgridError = HamiltonJacobi_ASGS_opt(coefficients,nd,stabFlag='2',lag=True)
    shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=rd_shockCapturingFactor,lag=True)
else:
    subgridError = HamiltonJacobi_ASGS_opt(coefficients,nd,stabFlag='2',lag=False)
    shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=rd_shockCapturingFactor,lag=False)

massLumping = False

numericalFluxType = DoNothing

multilevelNonlinearSolver  = MultilevelEikonalSolver
levelNonlinearSolver = UnstructuredFMMandFSWsolvers.FMMEikonalSolver
multilevelNonlinearSolver  = NLNI
levelNonlinearSolver = Newton

if rdtimeIntegration != 'newton':    
    maxLineSearches = 0
nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

#this needs to be set appropriately for pseudo-transient
tolFac = 0.95

if rdtimeIntegration != 'newton':
    nl_atol_res = 0.0
    maxNonlinearIts = 1
    maxLineSearches = 0
    psitc['nStepsForce']=5
    psitc['nStepsMax']=10
    nonlinearSolverConvergenceTest = 'rits'
    levelNonlinearSolverConvergenceTest = 'rits'
else:
    nl_atol_res = 0.0
    maxNonlinearIts = 10
    maxLineSearches = 0
    nonlinearSolverConvergenceTest = 'rits'
    levelNonlinearSolverConvergenceTest = 'rits'


matrix = SparseMatrix

if usePETSc:
    multilevelLinearSolver =PETSc # KSP_petsc4py
    levelLinearSolver = PETSc #KSP_petsc4py
    linear_solver_options_prefix = 'rdls_'
    linearSmoother = None
    #linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU
    
linTolFac = 0.001

conservativeFlux = None
