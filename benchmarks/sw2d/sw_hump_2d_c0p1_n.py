from proteus import *
from proteus.default_n import *
from sw_hump_2d_p import *

implicit=True

if implicit:
    timeIntegration = FLCBDF
    stepController = FLCBDF_controller
    timeIntegration = BackwardEuler
    stepController  = Min_dt_controller
    timeIntegration = BackwardEuler_cfl
    stepController = Min_dt_controller
    runCFL=0.33
    rtol_u[0] = 1.0e-4
    rtol_u[1] = 1.0e-4
    rtol_u[2] = 1.0e-4
    atol_u[0] = 1.0e-4
    atol_u[1] = 1.0e-4
    atol_u[2] = 1.0e-4
    femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
                 1:C0_AffineLinearOnSimplexWithNodalBasis,
                 2:C0_AffineLinearOnSimplexWithNodalBasis}
    elementQuadrature = SimplexGaussQuadrature(nd,3)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)
    #femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis,
    #             1:C0_AffineQuadraticOnSimplexWithNodalBasis,
    #             2:C0_AffineQuadraticOnSimplexWithNodalBasis}    
    # elementQuadrature = SimplexGaussQuadrature(nd,4)
    # elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)
    elementQuadrature = SimplexGaussQuadrature(nd,5)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,5)
    multilevelNonlinearSolver  = Newton
    
    levelNonlinearSolver = Newton

    fullNewtonFlag = True
    nDTout=201
else:
    runCFL=0.25
    timeOrder = 1
    class SSPRKwrap(LinearSSPRKintegration):
        def __init__(self,vt):
            LinearSSPRKintegration.__init__(self,vt,timeOrder,runCFL)
            return
    timeIntegration = SSPRKwrap 
    stepController=Min_dt_RKcontroller
    systemStepControllerType = SplitOperator.Sequential_MinAdaptiveModelStep
    #nDTout=101
    femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis,
                 1:DG_AffineP0_OnSimplexWithMonomialBasis,
                 2:DG_AffineP0_OnSimplexWithMonomialBasis}
#    femSpaces = {0:DG_Constants,
#                 1:DG_Constants,
#                 2:DG_Constants}
    elementQuadrature = SimplexGaussQuadrature(nd,1)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,1)
    numericalFluxType = ShallowWater_2D
    multilevelNonlinearSolver  = Newton
    levelNonlinearSolver = testStuff.SSPRKNewton
    fullNewtonFlag = True
#    femSpaces = {0:DG_AffineP1_OnSimplexWithMonomialBasis,
#                 1:DG_AffineP1_OnSimplexWithMonomialBasis}
#femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis,
#             1:C0_AffineQuadraticOnSimplexWithNodalBasis}

nnx=51
nny=51
he = L[0]/float(nnx-1)
triangleOptions="pAq30Dena%f"  % (0.5*he**2,)

#subgridError = ShallowWater_CFL(coefficients,nd,g)
subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,lag=False)
#added flag for using SUPG stabilization based on Berger and Stockstill, 95
try_supg_stabilization = False

massLumping=False

shockCapturing = None
shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False)

#numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG
numericalFluxType = ShallowWater_2D


tolFac = 0.0

nl_atol_res = 1.0e-4
l_atol_res = 1.0e-8
l_rtol_res = 0.0
#nl_atol_res = 1.0e-8

matrix = SparseMatrix
#matrix = numpy.array
multilevelLinearSolver = KSP_petsc4py

levelLinearSolver = KSP_petsc4py

linTolFac = 0.001

#conservativeFlux = {0:'pwl'}
