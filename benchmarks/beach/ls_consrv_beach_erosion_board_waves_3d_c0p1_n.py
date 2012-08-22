from proteus import *
from proteus.default_n import *
from beach_erosion_board_waves_3d import *
from ls_consrv_beach_erosion_board_waves_3d_p import *


timeIntegrator = ForwardIntegrator
timeIntegration = BackwardEuler#NoIntegration
timeIntegration = NoIntegration

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,sloshbox_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,sloshbox_quad_order)

subgridError = None
#subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=False)

massLumping = False

numericalFluxType = DoNothing

shockCapturing = None

#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False)
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False)

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-7#1.0e-10

maxNonlinearIts = 10

matrix = SparseMatrix


if usePETSc:
    multilevelLinearSolver = PETSc
    levelLinearSolver = PETSc
else:
    multilevelLinearSolver = LU

    levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 1.0e-6

conservativeFlux = None
