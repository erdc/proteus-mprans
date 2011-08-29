from proteus import *
from proteus.default_n import *
from ls_p import *

timeIntegration = BackwardEuler
stepController = FixedStep

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
elementQuadrature = SimplexGaussQuadrature(nd,quad_order)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order)

subgridError = HamiltonJacobi_ASGS_opt(coefficients,nd,lag=False)

massLumping = False

numericalFluxType =  Advection_DiagonalUpwind_IIPG_exterior

shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=ls_shockCapturingFactor,lag=False)

multilevelNonlinearSolver  = Newton
levelNonlinearSolver = Newton

nonlinearSmoother = None
linearSmoother = None

fullNewtonFlag = True

tolFac = 1e-3
nl_atol_res = 0.0

maxNonlinearIts = 10
maxLineSearches = 0

matrix = SparseMatrix

multilevelLinearSolver = PETSc
levelLinearSolver = PETSc


nonlinearSolverConvergenceTest = 'rits'
levelNonlinearSolverConvergenceTest = 'rits'
linear_solver_options_prefix = 'ncls_'

linTolFac = 0.001

conservativeFlux = None
