from proteus import *
from twp_navier_stokes_p import *
from wavetank import *

timeIntegration = BackwardEuler
stepController  = Min_dt_controller

femSpaces = {0:basis,
	     1:basis,
	     2:basis,
	     3:basis}

massLumping       = False
numericalFluxType = None
conservativeFlux  = None
numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior 
subgridError = NavierStokesASGS_velocity_pressure_optV2(coefficients,nd,lag=True,delayLagSteps=1,hFactor=hFactor,noPressureStabilization=False)
shockCapturing = NavierStokes_SC_opt(coefficients,nd,ns_shockCapturingFactor,lag=True)

fullNewtonFlag = True
multilevelNonlinearSolver = NewtonNS
levelNonlinearSolver      = NewtonNS

nonlinearSmoother = None
linearSmoother    = None

matrix = SparseMatrix

if useOldPETSc:
    multilevelLinearSolver = PETSc
    levelLinearSolver      = PETSc
else:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver      = KSP_petsc4py

if useSuperlu:
    multilevelLinearSolver = LU
    levelLinearSolver      = LU
    
linear_solver_options_prefix = 'rans2p_'
levelNonlinearSolverConvergenceTest = 'rits'
linearSolverConvergenceTest         = 'rits'

tolFac = 1e-3
nl_atol_res = 0.0

maxNonlinearIts = 10
maxLineSearches = 0
