#if True uses PETSc solvers
parallel = False
linearSmoother = None
#compute mass balance statistics or not
checkMass=False#True
#number of space dimensions
nd=2
#time integration, not relevant if using BDF with cfl timestepping
rtol_u = {0:1.0e-4}
atol_u = {0:1.0e-4}
rtol_res = {0:1.0e-4}
atol_res = {0:1.0e-4}
#
timeIntegration_vof = "BE"
timeIntegration_ls = "BE"
#if want bdf2 or bdf1
timeOrder = 1
runCFL = 0.3#0.3,0.185,0.125 for dgp1,dgp2,dgpk(3)
#
#spatial approximation orders
cDegree_ls=0 #0 -- CG. -1 -- DG
cDegree_vof=0
pDegree_ls=1 #level set 
pDegree_vof=pDegree_ls #volume of fluid should match ls for now
useHex=False#True
useMetrics=1.0
#
#spatial quadrature orders
#2*max(pDegree_vof,pDegree_ls)+1
if pDegree_ls == 2:
    vortex_quad_order = 5
else:
    vortex_quad_order = 3
#parallel partitioning info
from proteus import MeshTools
partitioningType = MeshTools.MeshParallelPartitioningTypes.node
#spatial mesh
lRefinement=3
#tag simulation name to level of refinement
#soname="vortexcgp2_bdf2_mc"+`lRefinement`
nn=nnx=nny=(2**lRefinement)*10+1
nnz=1
he=1.0/(nnx-1.0)
L=[1.0,1.0]

unstructured=False#True for tetgen, false for tet or hex from rectangular grid
if unstructured:
    from tank2dDomain import *
    domain = tank2d(L=L)
    bt = domain.boundaryTags
    domain.writePoly("tank2d")
    triangleOptions="pAq30Dena%8.8f"  % (0.5*he**2,)
else:
    from proteus.Domain import RectangularDomain
    domain = RectangularDomain(L)
#end time of simulation, full problem is T=8.0
T = 8.0#8.0#
#number of output time steps
nDTout = 80
#mass correction
applyCorrection=True
applyRedistancing=True
onlyVOF=False
#smoothing factors
#eps
epsFactHeaviside=1.5
epsFactDirac=1.5
epsFactDiffusion=10.0
epsFactRedistance=0.33
epsFact_vof=1.5
#
shockCapturingFactor_vof=0.1
shockCapturingFactor_ls=0.1
shockCapturingFactor_rd=0.1
lag_shockCapturing_vof=False
lag_shockCapturing_ls=False
lag_shockCapturing_rd=False
if useMetrics==0.0:
    shockCapturingFactor_vof=0.9
    shockCapturingFactor_ls=0.9
    shockCapturingFactor_rd=0.9
    lag_shockCapturing_vof=True
    lag_shockCapturing_ls=True
    lag_shockCapturing_rd=False

#use absolute tolerances on al models
atolRedistance = 0.1*he
atolConservation = 1.0e-6
atolVolumeOfFluid= 1.0e-6
atolLevelSet     = 1.0e-6
#controls 
linearSolverConvergenceTest = 'r-true' #rits is do a set number of iterations, r-true uses true residual, PETSc default is preconditioned residual
#redist solver
fmmFlag=0
#
#correctionType = 'dg'
#correctionType = 'dgp0'
#correctionType = 'global'
correctionType = 'cg'
#correctionType = 'none'
if useHex:
    hex=True
    soname="vortex_c0q"+`pDegree_ls`+correctionType+"_bdf_"+`timeOrder`+"_level_"+`lRefinement`
else:
    soname="vortex_c0p"+`pDegree_ls`+correctionType+"_bdf_"+`timeOrder`+"_level_"+`lRefinement`
