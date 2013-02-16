from proteus import *
from redist_p import *
from dtmb import *

#timeIntegration = BackwardEuler_cfl
#stepController = Osher_PsiTC_controller

#timeIntegration = BackwardEuler
#stepController = Osher_PsiTC_controller2	     

timeIntegrator  = ForwardIntegrator
timeIntegration = NoIntegration

femSpaces = {0:basis}
       
massLumping       = False
numericalFluxType = DoNothing    
conservativeFlux  = None
subgridError      = HamiltonJacobi_ASGS_opt(coefficients,nd,stabFlag='2',lag=False)
shockCapturing    = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=rd_shockCapturingFactor,lag=False)

fullNewtonFlag = True
multilevelNonlinearSolver  = NLNI
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
linearSolverConvergenceTest = 'r-true'

runCFL=1.0
rtol_res[0] = 0.0
atol_res[0] = 0.1*he
psitc['nStepsForce']=3
psitc['nStepsMax']=5
psitc['reduceRatio']=1.0
psitc['startRatio']=1.0 

tolFac = 0.0
nl_atol_res = 0.1*he
useEisenstatWalker = True

maxNonlinearIts = 5
maxLineSearches = 0
