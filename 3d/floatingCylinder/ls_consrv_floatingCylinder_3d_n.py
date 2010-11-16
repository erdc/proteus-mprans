from proteus import *
from proteus.default_n import *
from floatingCylinder import *
from ls_consrv_floatingCylinder_3d_p import *


stepController = FixedStep
timeIntegration = NoIntegration

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order)

subgridError = None

massLumping = False

numericalFluxType = DoNothing

if usePETSc:
    numericalFluxType = Diffusion_IIPG_exterior

shockCapturing = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 0.001*he#1.0e-10

maxNonlinearIts = 10

matrix = SparseMatrix

if usePETSc:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver = KSP_petsc4py
    linear_solver_options_prefix = 'mcorr_'
#    linearSmoother = StarILU
    linearSmoother = None
    linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU

linTolFac = 1.0e-6

conservativeFlux = None
