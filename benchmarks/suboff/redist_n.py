from proteus import *
from redist_p import *
from suboff import *

timeIntegration = BackwardEuler
stepController = Osher_PsiTC_controller2	     
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

if not use_petsc4py:
    multilevelLinearSolver = PETSc
    levelLinearSolver      = PETSc
else:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver      = KSP_petsc4py

linear_solver_options_prefix = 'rdls_'
nonlinearSolverConvergenceTest = 'rits'
linearSolverConvergenceTest = 'r-true'

runCFL=1.0
rtol_res[0] = 0.001
atol_res[0] = 0.0
psitc['nStepsForce']=0
psitc['nStepsMax']=0
psitc['reduceRatio']=1.0
psitc['startRatio']=1.0 

tolFac = 10.0
nl_atol_res = 0.0

maxNonlinearIts = 1
maxLineSearches = 0
