from proteus import *
from proteus.default_n import *
from ls_consrv_vortex_3d_p import *
from vortex import *


timeIntegrator = ForwardIntegrator
timeIntegration = NoIntegration
stepController = Newton_controller

if cDegree_ls==0:
    if useHex:
        if pDegree_ls==1:
            femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}
        elif pDegree_ls==2:
            femSpaces = {0:C0_AffineLagrangeOnCubeWithNodalBasis}
        elementQuadrature = CubeGaussQuadrature(nd,vortex_quad_order)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,vortex_quad_order)
    else:
        if pDegree_ls==1:
            femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
        elif pDegree_ls==2:
            femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}        
        elementQuadrature = SimplexGaussQuadrature(nd,vortex_quad_order)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,vortex_quad_order)
    if parallel or LevelModelType ==  MCorr.LevelModel:
        numericalFluxType = DoNothing#Diffusion_IIPG_exterior
elif cDegree_ls==-1:
    if pDegree_ls==0:
        femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis}
    elif pDegree_ls==1:
        femSpaces = {0:DG_AffineP1_OnSimplexWithMonomialBasis}
    elif pDegree_ls==2:
        femSpaces = {0:DG_AffineP2_OnSimplexWithMonomialBasis}
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_SIPG


# elementQuadrature = SimplexLobattoQuadrature(nd,1)

# elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

subgridError = None

massLumping = False

shockCapturing = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = atolConservation

maxNonlinearIts = 100

matrix = SparseMatrix

if parallel:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver = KSP_petsc4py
    linear_solver_options_prefix = 'mcorr_'
    linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU

    levelLinearSolver = LU

linTolFac = 1.0e-6

conservativeFlux = {}
if checkMass:
    auxiliaryVariables = [AuxiliaryVariables.ConservationHistoryMC("vortex3d"+`lRefinement`+"p"+`pDegree_ls`)]