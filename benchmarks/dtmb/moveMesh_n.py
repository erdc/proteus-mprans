from proteus import *
from proteus.default_n import *
from moveMesh_p import *

timeIntegration = NoIntegration

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,quad_order)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order)


nLevels = 1

subgridError = None

massLumping = False

numericalFluxType = Stress_IIPG_exterior

shockCapturing = None

multilevelNonlinearSolver  = Newton
levelNonlinearSolver = Newton

nonlinearSmoother = None
linearSmoother = None

fullNewtonFlag = True

tolFac = 1e-2

nl_atol_res = 0.0
maxNonlinearIts = 4#should be linear
maxLineSearches = 0

matrix = SparseMatrix

multilevelLinearSolver = PETSc
levelLinearSolver = PETSc
linear_solver_options_prefix = 'mesh_'
linearSmoother = None
linearSolverConvergenceTest = 'r-true'

linTolFac = 0.001

conservativeFlux = None

##auxiliaryVariables=[tro]
