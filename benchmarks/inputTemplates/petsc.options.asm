-rans2p_ksp_type bcgsl
-rans2p_ksp_atol 1.0e-6 -rans2p_ksp_rtol 0.0
-rans2p_ksp_monitor_true_residual 
-rans2p_pc_type asm -rans2p_pc_asm_type basic
-ncls_ksp_type   fgmres -ncls_pc_type   hypre -ncls_pc_hypre_type   boomeramg -ncls_ksp_atol  1.0e-6   -ncls_ksp_rtol  0.0 -ncls_ksp_monitor_true_residual
-vof_ksp_type    fgmres -vof_pc_type    hypre -vof_pc_hypre_type    boomeramg -vof_ksp_atol   1.0e-6   -vof_ksp_rtol   0.0 -vof_ksp_monitor_true_residual
-rdls_ksp_type   fgmres -rdls_pc_type   hypre -vof_pc_hypre_type    boomeramg -rdls_ksp_atol  1.0e-6   -rdls_ksp_rtol  0.0 -rdls_ksp_monitor_true_residual
-mcorr_ksp_type  cg     -mcorr_pc_type  hypre -mcorr_pc_hypre_type  boomeramg -mcorr_ksp_atol 1.0e-6   -mcorr_ksp_rtol 0.0 -mcorr_ksp_monitor_true_residual
-log_summary



