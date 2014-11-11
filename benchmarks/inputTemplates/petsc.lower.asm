#velocity is upper, pressure is lower
-rans2p_ksp_type fgmres
#-rans2p_ksp_atol 1.0e-4 -rans2p_ksp_rtol 0.0
#-rans2p_ksp_view
#-rans2p_ksp_monitor_true_residual
-rans2p_ksp_converged_reason
-rans2p_pc_type fieldsplit
# use this to discard the lower block of the LU factorization
#-rans2p_pc_fieldsplit_type schur -rans2p_pc_fieldsplit_schur_fact_type upper -rans2p_pc_fieldsplit_schur_precondition selfp
# use this to discard the upper block of the LU factorization
#-rans2p_pc_fieldsplit_type schur -rans2p_pc_fieldsplit_schur_fact_type lower -rans2p_pc_fieldsplit_schur_precondition selfp
# use full to do full Schur factorization (exact if velocity and pressure are solved exactly)

#-rans2p_pc_fieldsplit_type schur -rans2p_pc_fieldsplit_schur_fact_type full -rans2p_pc_fieldsplit_schur_precondition selfp

-rans2p_pc_fieldsplit_type schur -rans2p_pc_fieldsplit_schur_fact_type lower -rans2p_pc_fieldsplit_schur_precondition selfp

#-rans2p_fieldsplit_velocity_ksp_type fgmres
#-rans2p_fieldsplit_velocity_ksp_monitor
#-rans2p_fieldsplit_velocity_ksp_atol 1.0e-1 -rans2p_fieldsplit_velocity_ksp_rtol 1.0e-1
#-rans2p_fieldsplit_velocity_pc_type gamg
#-rans2p_fieldsplit_velocity_pc_type ilu


#asm solver for velocity block
#-rans2p_fieldsplit_velocity_ksp_monitor_true_residual
-rans2p_fieldsplit_velocity_ksp_type preonly#fgmres
#-rans2p_fieldsplit_velocity_ksp_atol 1.0e-10
#-rans2p_fieldsplit_velocity_ksp_rtol 0.0

#asm with direct on blocks for velocity block
-rans2p_fieldsplit_velocity_pc_type asm
-rans2p_fieldsplit_velocity_pc_asm_type basic
-rans2p_fieldsplit_velocity_sub_ksp_type preonly
-rans2p_fieldsplit_velocity_sub_pc_type lu
-rans2p_fieldsplit_velocity_sub_pc_factor_mat_solver_package superlu

#direct solver for  velocity block
#-rans2p_fieldsplit_velocity_ksp_type preonly
#-rans2p_fieldsplit_velocity_pc_type lu
#-rans2p_fieldsplit_velocity_pc_factor_mat_solver_package superlu_dist

#"direct" solver for pressure block
-rans2p_fieldsplit_pressure_ksp_type fgmres
#-rans2p_fieldsplit_pressure_ksp_monitor_true_residual
-rans2p_fieldsplit_pressure_ksp_atol 0.0
-rans2p_fieldsplit_pressure_ksp_rtol 0.01
-rans2p_fieldsplit_pressure_pc_type asm
-rans2p_fieldsplit_pressure_pc_asm_type basic
-rans2p_fieldsplit_pressure_sub_ksp_type preonly
-rans2p_fieldsplit_pressure_sub_pc_type lu
-rans2p_fieldsplit_pressure_sub_pc_factor_mat_solver_package superlu

#direct  solvers for the other models
-ncls_ksp_type   preonly -ncls_pc_type   lu -ncls_pc_factor_mat_solver_package   superlu_dist
-vof_ksp_type    preonly -vof_pc_type    lu -vof_pc_factor_mat_solver_package    superlu_dist
-rdls_ksp_type   preonly -rdls_pc_type   lu -rdls_pc_factor_mat_solver_package   superlu_dist
-mcorr_ksp_type  preonly -mcorr_pc_type  lu -mcorr_pc_factor_mat_solver_package  superlu_dist
-kappa_ksp_type preonly -kappa_pc_type lu -kappa_pc_factor_mat_solver_package superlu_dist
-dissipation_ksp_type preonly -dissipation_pc_type lu -dissipation_pc_factor_mat_solver_package superlu_dist
-log_summary