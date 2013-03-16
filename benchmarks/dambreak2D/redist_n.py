from proteus import *
from redist_p import *
from dambreak import *

#timeIntegration = BackwardEuler_cfl
#stepController = Osher_PsiTC_controller

#timeIntegration = BackwardEuler
#stepController = Osher_PsiTC_controller2	     

#timeIntegrator  = ForwardIntegrator
#timeIntegration = NoIntegration

timeIntegration = BackwardEuler_cfl
stepController = RDLS.PsiTC

runCFL=1.0
psitc['nStepsForce']=3
psitc['nStepsMax']=15
psitc['reduceRatio']=1.5
psitc['startRatio']=1.0

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
nonlinearSolverConvergenceTest = 'rits'
levelNonlinearSolverConvergenceTest = 'rits'
linearSolverConvergenceTest = 'r-true'

rtol_res[0] = 0.0
atol_res[0] = 0.1*he

tolFac = 0.01
nl_atol_res = 0.01*he
useEisenstatWalker = True

maxNonlinearIts = 1
maxLineSearches = 0
