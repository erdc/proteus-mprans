from proteus import *
from proteus.default_n import *
from wigley import *
from redist_wigley_3d_p import *

timeIntegration = BackwardEuler_cfl
stepController = Osher_PsiTC_controller2
runCFL=1.0
rtol_res[0] = 1e-3
atol_res[0] = 0.0

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order)

subgridError = HamiltonJacobi_ASGS_opt(coefficients,nd,stabFlag='2',lag=True)    
shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=rd_shockCapturingFactor,lag=True)

massLumping = False

numericalFluxType = DoNothing

multilevelNonlinearSolver  = MultilevelEikonalSolver
levelNonlinearSolver = UnstructuredFMMandFSWsolvers.FMMEikonalSolver
multilevelNonlinearSolver  = NLNI
levelNonlinearSolver = Newton

nl_atol_res = 10.0
maxNonlinearIts = 1
maxLineSearches = 0
psitc['nStepsForce']=5
psitc['nStepsMax']=10
nonlinearSolverConvergenceTest = 'rits'
levelNonlinearSolverConvergenceTest = 'rits'	
        
nonlinearSmoother = None
linearSmoother = None

matrix = SparseMatrix

multilevelLinearSolver = PETSc
levelLinearSolver = PETSc

linear_solver_options_prefix = 'rdls_'

linTolFac = 0.001

conservativeFlux = None
