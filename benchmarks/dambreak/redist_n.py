from proteus import *
from redist_p import *
from dambreak import *

if redist_Newton:
    timeIntegration = NoIntegration
    stepController = Newton_controller
    tolFac = 0.0
    nl_atol_res = 0.1*he
    useEisenstatWalker = True
    linTolFac = 0.0
    maxNonlinearIts = 10
    maxLineSearches = 0
    nonlinearSolverConvergenceTest = 'r'
    levelNonlinearSolverConvergenceTest = 'r'
    linearSolverConvergenceTest = 'r-true'
else:
    timeIntegration = BackwardEuler_cfl
    stepController = RDLS.PsiTC
    runCFL=0.33
    psitc['nStepsForce']=3
    psitc['nStepsMax']=15
    psitc['reduceRatio']=2.0
    psitc['startRatio']=1.0
    rtol_res[0] = 0.0
    atol_res[0] = 0.1*he
    tolFac = 0.0
    nl_atol_res = 0.01*he
    useEisenstatWalker = True
    linTolFac = 0.0
    maxNonlinearIts = 1
    maxLineSearches = 0
    nonlinearSolverConvergenceTest = 'rits'
    levelNonlinearSolverConvergenceTest = 'rits'
    linearSolverConvergenceTest = 'r-true'

femSpaces = {0:basis}
       
massLumping       = False
numericalFluxType = DoNothing    
conservativeFlux  = None
subgridError      = HamiltonJacobi_ASGS_opt(coefficients,nd,stabFlag='2',lag=False)
shockCapturing    = RDLS.ShockCapturing(coefficients,nd,shockCapturingFactor=rd_shockCapturingFactor,lag=False)

fullNewtonFlag = True
multilevelNonlinearSolver  = Newton
levelNonlinearSolver       = Newton

nonlinearSmoother = NLGaussSeidel
linearSmoother    = None

matrix = SparseMatrix

if useOldPETSc:
    multilevelLinearSolver = PETSc
    levelLinearSolver      = PETSc
else:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver      = KSP_petsc4py

if useSuperlu:
    multilevelLinearSolver = LU
    levelLinearSolver      = LU

linear_solver_options_prefix = 'rdls_'
