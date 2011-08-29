from proteus import *
from proteus.default_n import *
from DTMB import *
from vof_p import *

elementQuadrature = SimplexGaussQuadrature(nd,quad_order)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order)


timeIntegration = BackwardEuler
stepController  = FixedStep

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=vof_shockCapturingFactor,lag=True)#linear
subgridError = Advection_ASGS(coefficients=coefficients,nd=nd,lag=False)
massLumping = False
numericalFluxType = Advection_DiagonalUpwind_IIPG_exterior

multilevelNonlinearSolver  = Newton
levelNonlinearSolver = Newton

nonlinearSmoother = None
linearSmoother = None

fullNewtonFlag = True

tolFac = 1e-4

nl_atol_res = 0.0

maxNonlinearIts = 10
maxLineSearches = 0

matrix = SparseMatrix

multilevelLinearSolver = PETSc
levelLinearSolver = PETSc
linear_solver_options_prefix = 'vof_'

nonlinearSolverConvergenceTest = 'rits'
levelNonlinearSolverConvergenceTest = 'rits'

linTolFac = 0.001

conservativeFlux = None
