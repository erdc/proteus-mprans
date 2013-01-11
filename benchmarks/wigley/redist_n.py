from proteus import *
from proteus.default_n import *
from wigley import *
from redist_p import *

timeIntegration = BackwardEuler
stepController = Osher_PsiTC_controller2
runCFL=1.0
#rtol_res[0] = 1e-3
#atol_res[0] = 0.0

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order)

subgridError = HamiltonJacobi_ASGS_opt(coefficients,nd,stabFlag='2',lag=True)    
shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=False,lag=True)

massLumping = False

numericalFluxType = DoNothing

multilevelNonlinearSolver  = Newton
levelNonlinearSolver = Newton

nl_atol_res = 10.0
maxNonlinearIts = 1
maxLineSearches = 0
psitc['nStepsForce']=5
psitc['startRatio'] = 10.0 
psitc['reduceRatio']= 1.0
psitc['nStepsMax']=10
nonlinearSolverConvergenceTest = 'rits'
levelNonlinearSolverConvergenceTest = 'rits'	
        
nonlinearSmoother = None
linearSmoother = None

matrix = SparseMatrix

if not use_petsc4py:
    multilevelLinearSolver = PETSc
    levelLinearSolver      = PETSc
else:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver      = KSP_petsc4py

linear_solver_options_prefix = 'rdls_'

linTolFac = 0.001

conservativeFlux = None
