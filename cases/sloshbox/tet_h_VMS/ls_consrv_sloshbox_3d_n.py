from proteus import *
from proteus.default_n import *
from sloshbox3d import *
from ls_consrv_sloshbox_3d_p import *


timeIntegrator = ForwardIntegrator
timeIntegration = BackwardEuler#NoIntegration
timeIntegration = NoIntegration

if useHex:
    if spaceOrder == 1:
        femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}
    elif spaceOrder == 2:
        femSpaces = {0:C0_AffineLagrangeOnCubeWithNodalBasis}
    elementQuadrature = CubeGaussQuadrature(nd,sloshbox_quad_order)
    elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,sloshbox_quad_order)
else:
    if spaceOrder == 1:
        femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
    elif spaceOrder == 2:
        femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
    elementQuadrature = SimplexGaussQuadrature(nd,sloshbox_quad_order)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,sloshbox_quad_order)


subgridError = None

massLumping = False

numericalFluxType = DoNothing
shockCapturing = None

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

nonlinearSmoother = None

fullNewtonFlag = True

tolFac = 1e-4
nl_atol_res = 0.0

maxNonlinearIts = 10
maxLineSearches = 0

matrix = SparseMatrix


multilevelLinearSolver = PETSc
levelLinearSolver =  PETSc
linear_solver_options_prefix = 'mcorr_'
linearSmoother = None
linearSolverConvergenceTest = 'rits'


conservativeFlux = None
if checkMass:
    auxiliaryVariables = [AuxiliaryVariables.ConservationHistoryMC("sloshbox3d"+`lRefinement`+"p"+`spaceOrder`)]
