#ifndef RBLES_H
#define RBLES_H
#include <cmath>
#include <iostream>
#include "CompKernel.h"

#define nSpace 3
#define nQuadraturePoints_element 4
#define nDOF_trial_element 4
#define nDOF_mesh_trial_element 4
#define nDOF_test_element 4
#define nDOF_test_X_trial_element 16
#define nQuadraturePoints_elementBoundary 3
#define nElementBoundaries_element 4

namespace RBLES
{
  inline double smoothedHeaviside(double eps, double phi)
  {
    double H;
    if (phi > eps)
      H=1.0;
    else if (phi < -eps)
      H=0.0;
    else if (phi==0.0)
      H=0.5;
    else
      H = 0.5*(1.0 + phi/eps + sin(M_PI*phi/eps)/M_PI);
    return H;
  }

  inline double smoothedHeaviside_integral(double eps, double phi)
  {
    double HI;
    if (phi > eps)
      {
	HI= phi - eps							\
	  + 0.5*(eps + 0.5*eps*eps/eps - eps*cos(M_PI*eps/eps)/(M_PI*M_PI)) \
	  - 0.5*((-eps) + 0.5*(-eps)*(-eps)/eps - eps*cos(M_PI*(-eps)/eps)/(M_PI*M_PI));
      }
    else if (phi < -eps)
      {
	HI=0.0;
      }
    else
      {
	HI = 0.5*(phi + 0.5*phi*phi/eps - eps*cos(M_PI*phi/eps)/(M_PI*M_PI)) \
	  - 0.5*((-eps) + 0.5*(-eps)*(-eps)/eps - eps*cos(M_PI*(-eps)/eps)/(M_PI*M_PI));
      }
    return HI;
  }
 
  inline double smoothedDirac(double eps, double phi)
  {
    double d;
    if (phi > eps)
      d=0.0;
    else if (phi < -eps)
      d=0.0;
    else
      d = 0.5*(1.0 + cos(M_PI*phi/eps))/eps;
    return d;
  }


  inline
    void calculateSubgridError_tau(const double&  Ct_sge,
				   const double&  Cd_sge,
				   const double* G,
				   const double& G_dd_G,
				   const double& tr_G,
				   const double& rho,
				   const double& Dt,
				   const double v[nSpace],
				   const double& mu,
				   double& tau_p,
				   double& tau_v,
				   double& cfl)
  {
    const double rho2=rho*rho,Dt2=Dt*Dt,mu2=mu*mu;
    register double v_d_Gv=0.0;
    for(int I=0;I<nSpace;I++)
      for (int J=0;J<nSpace;J++)
	v_d_Gv += v[I]*G[I*nSpace+J]*v[J];
    cfl = 2.0*sqrt(v_d_Gv);
    //cek 1.0/sqrt(rho2*Dt2 + 4*v_d_Gv + 144*mu2*G_dd_G); ?
    /* tau_v = 1.0/sqrt(rho2*Dt2 + 4*v_d_Gv + 144*mu2*G_dd_G); */
    /* tau_p = 1.0/(tr_G*tau_v);  */
    //cek "correct" tau
    
    //std::cout<<rho2<<"  "<<Dt2<<"   "<< v_d_Gv<<std::endl;
    
    
    tau_v = 1.0/sqrt(Ct_sge*rho2*Dt2 + rho2*v_d_Gv + Cd_sge*mu2*G_dd_G);
    tau_p = 1.0/(tr_G*tau_v);
    //debug
    /* double tau_v_old = tau_v,tau_p_old = tau_p; */
    /* double nrm_v=0.0,h=1.0/20.0; */
    /* double oneByAbsdt =  fabs(Dt); */
    /* for(int I=0;I<nSpace;I++) */
    /* 	nrm_v += v[I]*v[I]; */
    /* nrm_v = sqrt(nrm_v); */
    /* cfl = nrm_v/h; */
    /* tau_v = 1.0/(4.0*mu/(h*h) + 2.0*rho*nrm_v/h + oneByAbsdt); */
    /* tau_p = 4.0*mu + 2.0*rho*nrm_v*h+ oneByAbsdt*h*h; */
    /* std::cout<<"nrm_v "<<nrm_v<<" tau_v "<<tau_v<<"\t"<<tau_v_old<<" tau_p "<<tau_p<<'\t'<<tau_p_old<<std::endl; */
  }

}//RBLES
extern "C"
{
  void calculateResidual_RBLES(//testing mesh replacement
				  double* mesh_trial_ref,
				  double* mesh_grad_trial_ref,
				  double* mesh_dof,
				  int* mesh_l2g,
				  double* dV_ref,
				  double* p_trial_ref,
				  double* p_grad_trial_ref,
				  double* p_test_ref,
				  double* p_grad_test_ref,
				  double* vel_trial_ref,
				  double* vel_grad_trial_ref,
				  double* vel_test_ref,
				  double* vel_grad_test_ref,
				  //element boundary
				  double* mesh_trial_trace_ref,
				  double* mesh_grad_trial_trace_ref,
				  double* dS_ref,
				  double* p_trial_trace_ref,
				  double* p_grad_trial_trace_ref,
				  double* p_test_trace_ref,
				  double* p_grad_test_trace_ref,
				  double* vel_trial_trace_ref,
				  double* vel_grad_trial_trace_ref,
				  double* vel_test_trace_ref,
				  double* vel_grad_test_trace_ref,					 
				  double* normal_ref,
				  double* boundaryJac_ref,
				  //end testing meshreplacement
				  int nElements_global,
                                  double useRBLES,
				  double alpha_bdf,
				  double eps_rho,
				  double eps_mu,
				  double sigma,
				  double rho_0,
				  double nu_0,
				  double rho_1,
				  double nu_1,
				  double Ct_sge,
				  double Cd_sge,
				  double C_dc,
				  int* p_l2g, int* vel_l2g,
				  double* p_dof, double* u_dof, double* v_dof, double* w_dof,
				  int* p_IBC, int* u_IBC, int* v_IBC, int* w_IBC,
				  double* g,
				  double* phi,
				  double* n,
				  double* kappa,
				  double* q_mom_u_acc,
				  double* q_mom_v_acc,
				  double* q_mom_w_acc,
				  double* q_mass_adv,
				  double* q_mom_u_acc_beta_bdf, double* q_mom_v_acc_beta_bdf, double* q_mom_w_acc_beta_bdf,
				  double* q_velocity_last,
				  double* q_cfl,
				  double* q_numDiff_u, double* q_numDiff_v, double* q_numDiff_w,
				  double* q_numDiff_u_last, double* q_numDiff_v_last, double* q_numDiff_w_last,
				  double* q_elementResidual_p, double* q_elementResidual_u, double* q_elementResidual_v, double* q_elementResidual_w,
				  int* sdInfo_u_u_rowptr,int* sdInfo_u_u_colind,			      
				  int* sdInfo_u_v_rowptr,int* sdInfo_u_v_colind,
				  int* sdInfo_u_w_rowptr,int* sdInfo_u_w_colind,
				  int* sdInfo_v_v_rowptr,int* sdInfo_v_v_colind,
				  int* sdInfo_v_u_rowptr,int* sdInfo_v_u_colind,
				  int* sdInfo_v_w_rowptr,int* sdInfo_v_w_colind,
				  int* sdInfo_w_w_rowptr,int* sdInfo_w_w_colind,
				  int* sdInfo_w_u_rowptr,int* sdInfo_w_u_colind,
				  int* sdInfo_w_v_rowptr,int* sdInfo_w_v_colind,
				  int offset_p, int offset_u, int offset_v, int offset_w, int stride_p, int stride_u, int stride_v, int stride_w, double* globalResidual,
				  int nExteriorElementBoundaries_global,
				  int* exteriorElementBoundariesArray,
				  int* elementBoundaryElementsArray,
				  int* elementBoundaryLocalElementBoundariesArray,
				  double* ebqe_phi_ext,
				  double* ebqe_n_ext,
				  double* ebqe_kappa_ext,
				  int* isDOFBoundary_p,
				  int* isDOFBoundary_u,
				  int* isDOFBoundary_v,
				  int* isDOFBoundary_w,
				  int* isAdvectiveFluxBoundary_p,
				  int* isAdvectiveFluxBoundary_u,
				  int* isAdvectiveFluxBoundary_v,
				  int* isAdvectiveFluxBoundary_w,
				  int* isDiffusiveFluxBoundary_u,
				  int* isDiffusiveFluxBoundary_v,
				  int* isDiffusiveFluxBoundary_w,
				  double* ebqe_bc_p_ext,
				  double* ebqe_bc_flux_mass_ext,
				  double* ebqe_bc_flux_mom_u_adv_ext,
				  double* ebqe_bc_flux_mom_v_adv_ext,
				  double* ebqe_bc_flux_mom_w_adv_ext,
				  double* ebqe_bc_u_ext,
				  double* ebqe_bc_flux_u_diff_ext,
				  double* ebqe_penalty_ext,
				  double* ebqe_bc_v_ext,
				  double* ebqe_bc_flux_v_diff_ext,
				  double* ebqe_bc_w_ext,
				  double* ebqe_bc_flux_w_diff_ext,
				  double* q_velocity,
				  double* ebqe_velocity_ext,
				  double* flux);

  void calculateJacobian_RBLES(//testing mesh replacement
				  double* mesh_trial_ref,
				  double* mesh_grad_trial_ref,
				  double* mesh_dof,
				  int* mesh_l2g,
				  double* dV_ref,
				  double* p_trial_ref,
				  double* p_grad_trial_ref,
				  double* p_test_ref,
				  double* p_grad_test_ref,
				  double* vel_trial_ref,
				  double* vel_grad_trial_ref,
				  double* vel_test_ref,
				  double* vel_grad_test_ref,
				  //element boundary
				  double* mesh_trial_trace_ref,
				  double* mesh_grad_trial_trace_ref,
				  double* dS_ref,
				  double* p_trial_trace_ref,
				  double* p_grad_trial_trace_ref,
				  double* p_test_trace_ref,
				  double* p_grad_test_trace_ref,
				  double* vel_trial_trace_ref,
				  double* vel_grad_trial_trace_ref,
				  double* vel_test_trace_ref,
				  double* vel_grad_test_trace_ref,					 
				  double* normal_ref,
				  double* boundaryJac_ref,
				  //end testing meshreplacement
				  int nElements_global,
                                  double useRBLES,
				  double alpha_bdf,
				  double eps_rho,
				  double eps_mu,
				  double sigma,
				  double rho_0,
				  double nu_0,
				  double rho_1,
				  double nu_1,
				  double Ct_sge,
				  double Cd_sge,
				  double C_dc,
				  int* p_l2g, int* vel_l2g,
				  double* p_dof, double* u_dof, double* v_dof, double* w_dof,
				  int* p_IBC, int* u_IBC, int* v_IBC, int* w_IBC,
				  double* g,
				  double* phi,
				  double* n,
				  double* kappa,
				  double* q_mom_u_acc_beta_bdf, double* q_mom_v_acc_beta_bdf, double* q_mom_w_acc_beta_bdf,
				  double* q_velocity_last,
				  double* q_cfl,
				  double* q_numDiff_u_last, double* q_numDiff_v_last, double* q_numDiff_w_last,
				  int* sdInfo_u_u_rowptr,int* sdInfo_u_u_colind,			      
				  int* sdInfo_u_v_rowptr,int* sdInfo_u_v_colind,
				  int* sdInfo_u_w_rowptr,int* sdInfo_u_w_colind,
				  int* sdInfo_v_v_rowptr,int* sdInfo_v_v_colind,
				  int* sdInfo_v_u_rowptr,int* sdInfo_v_u_colind,
				  int* sdInfo_v_w_rowptr,int* sdInfo_v_w_colind,
				  int* sdInfo_w_w_rowptr,int* sdInfo_w_w_colind,
				  int* sdInfo_w_u_rowptr,int* sdInfo_w_u_colind,
				  int* sdInfo_w_v_rowptr,int* sdInfo_w_v_colind,
				  int* csrRowIndeces_p_p,int* csrColumnOffsets_p_p,
				  int* csrRowIndeces_p_u,int* csrColumnOffsets_p_u,
				  int* csrRowIndeces_p_v,int* csrColumnOffsets_p_v,
				  int* csrRowIndeces_p_w,int* csrColumnOffsets_p_w,
				  int* csrRowIndeces_u_p,int* csrColumnOffsets_u_p,
				  int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
				  int* csrRowIndeces_u_v,int* csrColumnOffsets_u_v,
				  int* csrRowIndeces_u_w,int* csrColumnOffsets_u_w,
				  int* csrRowIndeces_v_p,int* csrColumnOffsets_v_p,
				  int* csrRowIndeces_v_u,int* csrColumnOffsets_v_u,
				  int* csrRowIndeces_v_v,int* csrColumnOffsets_v_v,
				  int* csrRowIndeces_v_w,int* csrColumnOffsets_v_w,
				  int* csrRowIndeces_w_p,int* csrColumnOffsets_w_p,
				  int* csrRowIndeces_w_u,int* csrColumnOffsets_w_u,
				  int* csrRowIndeces_w_v,int* csrColumnOffsets_w_v,
				  int* csrRowIndeces_w_w,int* csrColumnOffsets_w_w,
				  double* globalJacobian,
				  int nExteriorElementBoundaries_global,
				  int* exteriorElementBoundariesArray,
				  int* elementBoundaryElementsArray,
				  int* elementBoundaryLocalElementBoundariesArray,
				  double* ebqe_phi_ext,
				  double* ebqe_n_ext,
				  double* ebqe_kappa_ext,
				  int* isDOFBoundary_p,
				  int* isDOFBoundary_u,
				  int* isDOFBoundary_v,
				  int* isDOFBoundary_w,
				  int* isAdvectiveFluxBoundary_p,
				  int* isAdvectiveFluxBoundary_u,
				  int* isAdvectiveFluxBoundary_v,
				  int* isAdvectiveFluxBoundary_w,
				  int* isDiffusiveFluxBoundary_u,
				  int* isDiffusiveFluxBoundary_v,
				  int* isDiffusiveFluxBoundary_w,
				  double* ebqe_bc_p_ext,
				  double* ebqe_bc_flux_mass_ext,
				  double* ebqe_bc_flux_mom_u_adv_ext,
				  double* ebqe_bc_flux_mom_v_adv_ext,
				  double* ebqe_bc_flux_mom_w_adv_ext,
				  double* ebqe_bc_u_ext,
				  double* ebqe_bc_flux_u_diff_ext,
				  double* ebqe_penalty_ext,
				  double* ebqe_bc_v_ext,
				  double* ebqe_bc_flux_v_diff_ext,
				  double* ebqe_bc_w_ext,
				  double* ebqe_bc_flux_w_diff_ext,
				  int* csrColumnOffsets_eb_p_p,
				  int* csrColumnOffsets_eb_p_u,
				  int* csrColumnOffsets_eb_p_v,
				  int* csrColumnOffsets_eb_p_w,
				  int* csrColumnOffsets_eb_u_p,
				  int* csrColumnOffsets_eb_u_u,
				  int* csrColumnOffsets_eb_u_v,
				  int* csrColumnOffsets_eb_u_w,
				  int* csrColumnOffsets_eb_v_p,
				  int* csrColumnOffsets_eb_v_u,
				  int* csrColumnOffsets_eb_v_v,
				  int* csrColumnOffsets_eb_v_w,
				  int* csrColumnOffsets_eb_w_p,
				  int* csrColumnOffsets_eb_w_u,
				  int* csrColumnOffsets_eb_w_v,
				  int* csrColumnOffsets_eb_w_w);
  void calculateVelocityAverage_RBLES(int* permutations,
					 int nExteriorElementBoundaries_global,
					 int* exteriorElementBoundariesArray,
					 int nInteriorElementBoundaries_global,
					 int* interiorElementBoundariesArray,
					 int* elementBoundaryElementsArray,
					 int* elementBoundaryLocalElementBoundariesArray,
					 int* vel_l2g, 
					 double* u_dof, double* v_dof, double* w_dof,
					 double* vel_trial,
					 double* ebqe_velocity,
					 double* velocityAverage);
}//extern "C"
#endif
