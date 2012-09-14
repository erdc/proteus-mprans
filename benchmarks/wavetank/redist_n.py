from proteus import *
from redist_p import *
from wavetank import *

timeIntegration = BackwardEuler
stepController = Osher_PsiTC_controller2	     

timeIntegration = NoIntegration
stepController  = Newton_controller
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
linearSolverConvergenceTest = 'rits'

runCFL=1.0
rtol_res[0] = 0.0
atol_res[0] = 0.1*he
psitc['nStepsForce']=3
psitc['nStepsMax']=10 
psitc['reduceRatio']=0.5
psitc['startRatio']=1.0 

tolFac = 0.0
nl_atol_res = 0.1*he

maxNonlinearIts = 1
maxLineSearches = 0
