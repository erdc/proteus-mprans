from proteus import *
from proteus.default_n import *
from redist_obstacleInTank_3d_p import *
from obstacleInTank3d import *

if rdtimeIntegration == 'newton':    
    timeIntegration = NoIntegration
    stepController = Newton_controller
else:
    timeIntegration = BackwardEuler_cfl
    stepController = Osher_PsiTC_controller
    runCFL=1.0
    rtol_res[0] = 0.0
    atol_res[0] = he*0.01#1.0e-6
    psitc['nStepsForce']=3
    psitc['nStepsMax']=5
if spaceOrder == 1:
    femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
if spaceOrder == 2:
    femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,obstacleInTank_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,obstacleInTank_quad_order)

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
tolFac = 0.0

if rdtimeIntegration != 'newton':
    nl_atol_res = 0.1*he
    maxNonlinearIts = 1
    maxLineSearches =0
    nonlinearSolverConvergenceTest = 'its'
    levelNonlinearSolverConvergenceTest = 'its'
else:
    nl_atol_res = 0.01*he#1.0e-7#0.01*L[0]/nnx
    maxNonlinearIts = 50
    maxLineSearches =100


matrix = SparseMatrix

if usePETSc:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver = KSP_petsc4py
    #multilevelLinearSolver = PETSc
    #levelLinearSolver = PETSc
    linear_solver_options_prefix = 'rdls_'
    linearSmoother = StarILU
#    linearSmoother = None
    linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU
    
linTolFac = 0.001

conservativeFlux = None
