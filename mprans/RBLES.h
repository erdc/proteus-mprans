#ifndef RBLES_H
#define RBLES_H
#include <cmath>
#include <iostream>
#include "CompKernel.h"
#include "ModelFactory.h"
//#define COMPRESSIBLE_FORM
namespace proteus
{
  class RBLES_base
  {
  public:
    virtual ~RBLES_base(){}
    virtual void calculateResidual(//element
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
				   //physics
				   double* elementDiameter,
				   double hFactor,
				   int nElements_global,
				   double useRBLES,
			           double useMetrics, 
				   double alphaBDF,
				   double epsFact_rho,
				   double epsFact_mu, 
				   double sigma,
				   double rho_0,
				   double nu_0,
				   double rho_1,
				   double nu_1,
				   double Ct_sge,
				   double Cd_sge,
				   double C_dc,
				   int* p_l2g, 
				   int* vel_l2g, 
				   double* p_dof, 
				   double* u_dof, 
				   double* v_dof, 
				   double* w_dof,
				   double* g,
				   double* phi,
				   double* normal_phi,
				   double* kappa_phi,
				   double* q_mom_u_acc,
				   double* q_mom_v_acc,
				   double* q_mom_w_acc,
				   double* q_mass_adv,
				   double* q_mom_u_acc_beta_bdf, double* q_mom_v_acc_beta_bdf, double* q_mom_w_acc_beta_bdf,
				   double* q_velocity_sge,
				   double* q_cfl,
				   double* q_numDiff_u, double* q_numDiff_v, double* q_numDiff_w,
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
				   int offset_p, int offset_u, int offset_v, int offset_w, 
				   int stride_p, int stride_u, int stride_v, int stride_w, 
				   double* globalResidual,
				   int nExteriorElementBoundaries_global,
				   int* exteriorElementBoundariesArray,
				   int* elementBoundaryElementsArray,
				   int* elementBoundaryLocalElementBoundariesArray,
				   double* ebqe_phi_ext,
				   double* ebqe_normal_phi_ext,
				   double* ebqe_kappa_phi_ext,
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
				   double* ebqe_velocity,
				   double* flux,
				   double* elementResidual_p)=0;
    virtual void calculateJacobian(//element
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
				   //physics
				   double* elementDiameter,
				   double hFactor,
				   int nElements_global,
				   double useRBLES,
			           double useMetrics, 
				   double alphaBDF,
				   double epsFact_rho,
				   double epsFact_mu,
				   double sigma,
				   double rho_0,
				   double nu_0,
				   double rho_1,
				   double nu_1,
				   double Ct_sge,
				   double Cd_sge,
				   double C_dg,
				   int* p_l2g, 
				   int* vel_l2g,
				   double* p_dof, double* u_dof, double* v_dof, double* w_dof,
				   double* g,
				   double* phi,
				   double* normal_phi,
				   double* kappa_phi,
				   double* q_mom_u_acc_beta_bdf, double* q_mom_v_acc_beta_bdf, double* q_mom_w_acc_beta_bdf,
				   double* q_velocity_sge,
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
				   double* ebqe_normal_phi_ext,
				   double* ebqe_kappa_phi_ext,
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
				   int* csrColumnOffsets_eb_w_w)=0;
    virtual void calculateVelocityAverage(int nExteriorElementBoundaries_global,
    					  int* exteriorElementBoundariesArray,
    					  int nInteriorElementBoundaries_global,
    					  int* interiorElementBoundariesArray,
    					  int* elementBoundaryElementsArray,
    					  int* elementBoundaryLocalElementBoundariesArray,
    					  double* mesh_dof,
    					  int* mesh_l2g,
    					  double* mesh_trial_trace_ref,
    					  double* mesh_grad_trial_trace_ref,
    					  double* normal_ref,
    					  double* boundaryJac_ref,
    					  int* vel_l2g,
    					  double* u_dof,
    					  double* v_dof,
    					  double* w_dof,
    					  double* vel_trial_trace_ref,
    					  double* ebqe_velocity,
    					  double* velocityAverage)=0;
  };
  
  template<class CompKernelType,
	   int nSpace,
	   int nQuadraturePoints_element,
	   int nDOF_mesh_trial_element,
	   int nDOF_trial_element,
	   int nDOF_test_element,
	   int nQuadraturePoints_elementBoundary>
  class RBLES : public RBLES_base
  {
  public:
    const int nDOF_test_X_trial_element;
    CompKernelType ck;
    RBLES():
      nDOF_test_X_trial_element(nDOF_test_element*nDOF_trial_element),
      ck()
    {}

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
    void calculateSubgridError_tau(const double&  hFactor,
				   const double& elementDiameter,
				   const double& dmt,
				   const double& dm,
				   const double df[nSpace],
				   const double& a,
				   double& tau_v,
				   double& tau_p,
				   double& cfl)
    {
      double h,oneByAbsdt,density,viscosity,nrm_df;
      h = hFactor*elementDiameter;
      density = dm;
      viscosity =  a;
      nrm_df=0.0;
      for(int I=0;I<nSpace;I++)
	nrm_df+=df[I]*df[I];
      nrm_df = sqrt(nrm_df);
      cfl = nrm_df/(h*density);
      oneByAbsdt =  fabs(dmt);
      tau_v = 1.0/(4.0*viscosity/(h*h) + 2.0*nrm_df/h + oneByAbsdt);
      tau_p = 4.0*viscosity + 2.0*nrm_df*h + oneByAbsdt*h*h;
    }


    inline
    void calculateSubgridError_tau(     const double&  Ct_sge,
                                        const double&  Cd_sge,
			                const double   G[nSpace*nSpace],
					const double&  G_dd_G,
					const double&  tr_G,
					const double&  A0,
					const double   Ai[nSpace],
					const double&  Kij,
					double& tau_v,
					double& tau_p,
					double& q_cfl)	
    {
      double v_d_Gv=0.0; 
      for(int I=0;I<nSpace;I++) 
         for (int J=0;J<nSpace;J++) 
           v_d_Gv += Ai[I]*G[I*nSpace+J]*Ai[J];     
    
      tau_v = 1.0/sqrt(Ct_sge*A0*A0 + v_d_Gv + Cd_sge*Kij*Kij*G_dd_G); 
      tau_p = 1.0/(tr_G*tau_v);     
    }


    void calculateResidual(//element
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
			   //physics
			   double* elementDiameter,
			   double hFactor,
			   int nElements_global,
			   double useRBLES,
			   double useMetrics, 
			   double alphaBDF,
			   double epsFact_rho,
			   double epsFact_mu, 
			   double sigma,
			   double rho_0,
			   double nu_0,
			   double rho_1,
			   double nu_1,
			   double Ct_sge,
			   double Cd_sge,
			   double C_dc,
			   int* p_l2g, 
			   int* vel_l2g, 
			   double* p_dof, 
			   double* u_dof, 
			   double* v_dof, 
			   double* w_dof,
			   double* g,
			   double* phi,
			   double* normal_phi,
			   double* kappa_phi,
			   double* q_mom_u_acc,
			   double* q_mom_v_acc,
			   double* q_mom_w_acc,
			   double* q_mass_adv,
			   double* q_mom_u_acc_beta_bdf, double* q_mom_v_acc_beta_bdf, double* q_mom_w_acc_beta_bdf,
			   double* q_velocity_sge,
			   double* q_cfl,
			   double* q_numDiff_u, double* q_numDiff_v, double* q_numDiff_w,
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
			   int offset_p, int offset_u, int offset_v, int offset_w, 
			   int stride_p, int stride_u, int stride_v, int stride_w, 
			   double* globalResidual,
			   int nExteriorElementBoundaries_global,
			   int* exteriorElementBoundariesArray,
			   int* elementBoundaryElementsArray,
			   int* elementBoundaryLocalElementBoundariesArray,
			   double* ebqe_phi_ext,
			   double* ebqe_normal_phi_ext,
			   double* ebqe_kappa_phi_ext,
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
			   double* ebqe_velocity,
			   double* flux,
			   double* elementResidual_p_save)
    {
  //double globalConservationError=0.0;
  for(int eN=0;eN<nElements_global;eN++)
    {
      //declare local storage for element residual and initialize
      register double elementResidual_p[nDOF_test_element],
	elementResidual_u[nDOF_test_element],
	elementResidual_v[nDOF_test_element],
	elementResidual_w[nDOF_test_element],
	eps_rho,eps_mu;
      for (int i=0;i<nDOF_test_element;i++)
	{
	  elementResidual_p[i]=0.0;
	  elementResidual_u[i]=0.0;
	  elementResidual_v[i]=0.0;
	  elementResidual_w[i]=0.0;
	}//i
      //
      //loop over quadrature points and compute integrands
      //
      for(int k=0;k<nQuadraturePoints_element;k++)
        {
	  //compute indices and declare local storage
	  register int eN_k = eN*nQuadraturePoints_element+k,
	    eN_k_nSpace = eN_k*nSpace,
	    eN_nDOF_trial_element = eN*nDOF_trial_element;
	  register double p=0.0,u=0.0,v=0.0,w=0.0,
	    grad_p[nSpace],grad_u[nSpace],grad_v[nSpace],grad_w[nSpace],
	    mom_u_acc=0.0,
	    dmom_u_acc_u=0.0,
	    mom_v_acc=0.0,
	    dmom_v_acc_v=0.0,
	    mom_w_acc=0.0,
	    dmom_w_acc_w=0.0,
	    mom_u_acc_t=0.0,
	    dmom_u_acc_u_t=0.0,
	    mom_v_acc_t=0.0,
	    dmom_v_acc_v_t=0.0,
	    mom_w_acc_t=0.0,
	    dmom_w_acc_w_t=0.0,
	    pdeResidual_p=0.0,
	    pdeResidual_u=0.0,
	    pdeResidual_v=0.0,
	    pdeResidual_w=0.0,
	    subgrid_p=0.0,
	    subgrid_u=0.0,
	    subgrid_v=0.0,
	    subgrid_w=0.0, 
	    tau_0=0.0,
	    tau_1=0.0,
	    jac[nSpace*nSpace],
	    jacDet,
	    jacInv[nSpace*nSpace],
	    p_grad_trial[nDOF_trial_element*nSpace],vel_grad_trial[nDOF_trial_element*nSpace],
	    p_test_dV[nDOF_trial_element],vel_test_dV[nDOF_trial_element],
	    p_grad_test_dV[nDOF_test_element*nSpace],vel_grad_test_dV[nDOF_test_element*nSpace],
	    dV,x,y,z,
	    G[nSpace*nSpace],G_dd_G,tr_G,norm_Rv,h_phi,
	    rho,mu,H_rho,d_rho,H_mu,d_mu;
	  //get jacobian, etc for mapping reference element
	  ck.calculateMapping_element(eN,
				      k,
				      mesh_dof,
				      mesh_l2g,
				      mesh_trial_ref, 
				      mesh_grad_trial_ref,
				      jac,
				      jacDet,
				      jacInv,
				      x,y,z);
	  //get the physical integration weight
	  dV = fabs(jacDet)*dV_ref[k];
	  ck.calculateG(jacInv,G,G_dd_G,tr_G);
	  ck.calculateGScale(G,&normal_phi[eN_k_nSpace],h_phi);

	  eps_rho = epsFact_rho*h_phi;
	  eps_mu = epsFact_mu*h_phi;
	  
	      eps_rho = epsFact_rho*elementDiameter[eN];
	      eps_mu  = epsFact_mu *elementDiameter[eN];
	  
	  
	  //get the trial function gradients
	  ck.gradTrialFromRef(&p_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,p_grad_trial);
	  ck.gradTrialFromRef(&vel_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,vel_grad_trial);
	  //get the solution
	  ck.valFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],&p_trial_ref[k*nDOF_trial_element],p);
	  ck.valFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],u);
	  ck.valFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],v);
	  ck.valFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],w);
	  //get the solution gradients
	  ck.gradFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],p_grad_trial,grad_p);
	  ck.gradFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_u);
	  ck.gradFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_v);
	  ck.gradFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_w);
	  //precalculate test function products with integration weights
	  for (int j=0;j<nDOF_trial_element;j++)
	    {
	      p_test_dV[j] = p_test_ref[k*nDOF_trial_element+j]*dV;
	      vel_test_dV[j] = vel_test_ref[k*nDOF_trial_element+j]*dV;
	      for (int I=0;I<nSpace;I++)
		{
		  p_grad_test_dV[j*nSpace+I]   = p_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
		  vel_grad_test_dV[j*nSpace+I] = vel_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
		}
	    }
	  
		  
          H_rho = RBLES::smoothedHeaviside(eps_rho,phi[eN_k]);
          d_rho = RBLES::smoothedDirac(eps_rho,phi[eN_k]);
          H_mu  = RBLES::smoothedHeaviside(eps_mu,phi[eN_k]);
          d_mu  = RBLES::smoothedDirac(eps_mu,phi[eN_k]);
  
          rho = rho_0*(1.0-H_rho)+rho_1*H_rho;
          mu  = nu_0*(1.0-H_mu)+nu_1*H_mu;
 
           //u momentum accumulation
          mom_u_acc=u;
          dmom_u_acc_u=1.0;
  
          //v momentum accumulation
          mom_v_acc=v;
          dmom_v_acc_v=1.0;
  
          //w momentum accumulation
          mom_w_acc=w;
          dmom_w_acc_w=1.0;
          
	  //
	  //save momentum for time history 
	  //
	  q_mom_u_acc[eN_k] = mom_u_acc;			    
	  q_mom_v_acc[eN_k] = mom_v_acc;			    
	  q_mom_w_acc[eN_k] = mom_w_acc;

          // Interpolate velocity for passing down
          q_velocity[eN_k_nSpace+0]= u;
          q_velocity[eN_k_nSpace+1]= v;
          q_velocity[eN_k_nSpace+2]= w;  

          //
          //calculate time derivative at quadrature points
          //
          ck.bdf(alphaBDF,
		 q_mom_u_acc_beta_bdf[eN_k],
		 mom_u_acc,
		 dmom_u_acc_u,
		 mom_u_acc_t,
		 dmom_u_acc_u_t);
          ck.bdf(alphaBDF,
		 q_mom_v_acc_beta_bdf[eN_k],
		 mom_v_acc,
		 dmom_v_acc_v,
		 mom_v_acc_t,
		 dmom_v_acc_v_t);
          ck.bdf(alphaBDF,
		 q_mom_w_acc_beta_bdf[eN_k],
		 mom_w_acc,
		 dmom_w_acc_w,
		 mom_w_acc_t,
		 dmom_w_acc_w_t);

          double Dt=(dmom_u_acc_u_t/dmom_u_acc_u);
          double v_d_Gv= 0.0;
	  double vel[nSpace] = {u,v,w};
          for(int I=0;I<nSpace;I++)
            for (int J=0;J<nSpace;J++)
	      v_d_Gv += vel[I]*G[I*nSpace+J]*vel[J];
  
          tau_0 = 1.0/sqrt(Ct_sge*rho*rho*Dt*Dt + rho*rho*v_d_Gv + Cd_sge*mu*mu*G_dd_G);
          tau_0 = 1.0/sqrt(Ct_sge*rho*rho*Dt*Dt  + Cd_sge*mu*mu*G_dd_G);
          tau_1 = 1.0/(tr_G*tau_0);
          q_cfl[eN_k] = 2.0*sqrt(v_d_Gv);						

          //
          //calculate strong residual
          //
	  pdeResidual_p = grad_u[0] + grad_v[1] + grad_w[2];

	  pdeResidual_u = rho *(mom_u_acc_t + u*grad_u[0] + v*grad_u[1] + w*grad_u[2] - g[0])
	                + grad_p[0]; // - mu poisson_u  ==> ommited as only first order right now !! 

	  pdeResidual_v = rho *(mom_v_acc_t + u*grad_v[0] + v*grad_v[1] + w*grad_v[2] - g[1])
	                + grad_p[1]; // - mu poisson_u  ==> ommited as only first order right now !! 
	                
	  pdeResidual_w = rho *(mom_w_acc_t + u*grad_w[0] + v*grad_w[1] + w*grad_w[2] - g[2])
	                + grad_p[2]; // - mu poisson_u  ==> ommited as only first order right now !! 

	  pdeResidual_u = rho *(mom_u_acc_t - g[0])
	                + grad_p[0]; // - mu poisson_u  ==> ommited as only first order right now !! 

	  pdeResidual_v = rho *(mom_v_acc_t - g[1])
	                + grad_p[1]; // - mu poisson_u  ==> ommited as only first order right now !! 
	                
	  pdeResidual_w = rho *(mom_w_acc_t - g[2])
	                + grad_p[2]; // - mu poisson_u  ==> ommited as only first order right now !! 

          //
          //calculate reconstructed subgrid scales
          //
          subgrid_u = -tau_0*pdeResidual_u;
          subgrid_v = -tau_0*pdeResidual_v;
          subgrid_w = -tau_0*pdeResidual_w;

          subgrid_p = -tau_1*pdeResidual_p;

          //
          //calculate discontinuity capturing addition to viscosity
          //
	  norm_Rv = sqrt(pdeResidual_u*pdeResidual_u + pdeResidual_v*pdeResidual_v + pdeResidual_w*pdeResidual_w);
	  //mu += C_dc*norm_Rv/sqrt(G_dd_G);

          // 
          //update element residual 
          // 
          for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      register int i_nSpace=i*nSpace;	
				    

	      elementResidual_p[i] += -pdeResidual_p*p_test_dV[i]
	                            + subgrid_u*p_grad_test_dV[i_nSpace+0]
	                            + subgrid_v*p_grad_test_dV[i_nSpace+1]
	                            + subgrid_w*p_grad_test_dV[i_nSpace+2];

	      elementResidual_u[i] += rho *(mom_u_acc_t - g[0])*vel_test_dV[i]	//rho *(mom_u_acc_t+ u*grad_u[0] + v*grad_u[1] + w*grad_u[2] - g[0])*vel_test_dV[i]				        
			               + mu*(grad_u[0] + grad_u[0])*vel_grad_test_dV[i_nSpace+0]
			               + mu*(grad_u[1] + grad_v[0])*vel_grad_test_dV[i_nSpace+1]				   
			               + mu*(grad_u[2] + grad_w[0])*vel_grad_test_dV[i_nSpace+2]	
	                               -(p+subgrid_p)*vel_grad_test_dV[i_nSpace+0]     	 
		                    ;//    -rho*vel_grad_test_dV[i_nSpace+0] * (subgrid_u*u) //+ u*subgrid_u + subgrid_u*subgrid_u)
		                     //   -rho*vel_grad_test_dV[i_nSpace+1] * (subgrid_u*v) //+ u*subgrid_v + subgrid_u*subgrid_v)
		                     //   -rho*vel_grad_test_dV[i_nSpace+2] * (subgrid_u*w);// + u*subgrid_w + subgrid_u*subgrid_w);
		            
	      elementResidual_v[i] +=  rho *(mom_v_acc_t - g[1])*vel_test_dV[i]	//rho *(mom_v_acc_t+ u*grad_v[0] + v*grad_v[1] + w*grad_v[2] - g[1])*vel_test_dV[i]
			               + mu*(grad_v[0] + grad_u[1])*vel_grad_test_dV[i_nSpace+0]
			               + mu*(grad_v[1] + grad_v[1])*vel_grad_test_dV[i_nSpace+1]				   
			               + mu*(grad_v[2] + grad_w[1])*vel_grad_test_dV[i_nSpace+2]	
	                               -(p+subgrid_p)*vel_grad_test_dV[i_nSpace+1] 	                                  
		                    ;//   -rho*vel_grad_test_dV[i_nSpace+0] * (subgrid_v*u) //+ v*subgrid_u + subgrid_v*subgrid_u)
		                     //  -rho*vel_grad_test_dV[i_nSpace+1] * (subgrid_v*v) //+ v*subgrid_v + subgrid_v*subgrid_v)
		                     //  -rho*vel_grad_test_dV[i_nSpace+2] * (subgrid_v*w);// + v*subgrid_w + subgrid_v*subgrid_w);


	      elementResidual_w[i] +=  rho *(mom_w_acc_t - g[2])*vel_test_dV[i]	//rho *(mom_w_acc_t + u*grad_w[0] + v*grad_w[1] + w*grad_w[2] - g[2])*vel_test_dV[i]
			               + mu*(grad_w[0] + grad_u[2])*vel_grad_test_dV[i_nSpace+0]
			               + mu*(grad_w[1] + grad_v[2])*vel_grad_test_dV[i_nSpace+1]				   
			               + mu*(grad_w[2] + grad_w[2])*vel_grad_test_dV[i_nSpace+2]		 
	                               -(p+subgrid_p)*vel_grad_test_dV[i_nSpace+2] 
		                   ; //    -rho*vel_grad_test_dV[i_nSpace+0] * (subgrid_w*u) //+ w*subgrid_u + subgrid_w*subgrid_u)
		                    //    -rho*vel_grad_test_dV[i_nSpace+1] * (subgrid_w*v) //+ w*subgrid_v + subgrid_w*subgrid_v)
		                    //    -rho*vel_grad_test_dV[i_nSpace+2] * (subgrid_w*w);// + w*subgrid_w + subgrid_w*subgrid_w);
           }//i
	}
	  //
	  //load element into global residual and save element residual
	  //
	  for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      register int eN_i=eN*nDOF_test_element+i;

	      elementResidual_p_save[eN_i] +=  elementResidual_p[i];
	  
	      globalResidual[offset_p+stride_p*p_l2g[eN_i]]+=elementResidual_p[i];
	      globalResidual[offset_u+stride_u*vel_l2g[eN_i]]+=elementResidual_u[i];
	      globalResidual[offset_v+stride_v*vel_l2g[eN_i]]+=elementResidual_v[i];
	      globalResidual[offset_w+stride_w*vel_l2g[eN_i]]+=elementResidual_w[i];
	    }//i
	}//elements
      //
      //loop over exterior element boundaries to calculate surface integrals and load into element and global residuals
      //
      //ebNE is the Exterior element boundary INdex
      //ebN is the element boundary INdex
      //eN is the element index
     for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
	{ 
	  register int ebN = exteriorElementBoundariesArray[ebNE], 
	    eN  = elementBoundaryElementsArray[ebN*2+0],
	    ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0],
	    eN_nDOF_trial_element = eN*nDOF_trial_element;
	  register double elementResidual_p[nDOF_test_element],
	    elementResidual_u[nDOF_test_element],
	    elementResidual_v[nDOF_test_element],
	    elementResidual_w[nDOF_test_element],
	    eps_rho,eps_mu;
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      elementResidual_p[i]=0.0;
	      elementResidual_u[i]=0.0;
	      elementResidual_v[i]=0.0;
	      elementResidual_w[i]=0.0;
	    }
	  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	    { 
	      register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
		ebNE_kb_nSpace = ebNE_kb*nSpace,
		ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
		ebN_local_kb_nSpace = ebN_local_kb*nSpace;
	      register double p_ext=0.0,
		u_ext=0.0,
		v_ext=0.0,
		w_ext=0.0,
		grad_p_ext[nSpace],
		grad_u_ext[nSpace],
		grad_v_ext[nSpace],
		grad_w_ext[nSpace],

		jac_ext[nSpace*nSpace],
		jacDet_ext,
		jacInv_ext[nSpace*nSpace],
		boundaryJac[nSpace*(nSpace-1)],
		metricTensor[(nSpace-1)*(nSpace-1)],
		metricTensorDetSqrt,
		dS,p_test_dS[nDOF_test_element],vel_test_dS[nDOF_test_element],vel_grad_test_dS[nDOF_test_element*nSpace],
		p_grad_trial_trace[nDOF_trial_element*nSpace],vel_grad_trial_trace[nDOF_trial_element*nSpace],
		normal[3],x_ext,y_ext,z_ext,
		G[nSpace*nSpace],G_dd_G,tr_G,h_phi,h_penalty,penalty,gi[nSpace], unormal,gnormal,uneg;
	      //compute information about mapping from reference element to physical element
	      ck.calculateMapping_elementBoundary(eN,
						  ebN_local,
						  kb,
						  ebN_local_kb,
						  mesh_dof,
						  mesh_l2g,
						  mesh_trial_trace_ref,
						  mesh_grad_trial_trace_ref,
						  boundaryJac_ref,
						  jac_ext,
						  jacDet_ext,
						  jacInv_ext,
						  boundaryJac,
						  metricTensor,
						  metricTensorDetSqrt,
						  normal_ref,
						  normal,
						  x_ext,y_ext,z_ext);
	      dS = metricTensorDetSqrt*dS_ref[kb];
	      //get the metric tensor
	      //cek todo use symmetry
	      ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
	      ck.calculateGScale(G,&ebqe_normal_phi_ext[ebNE_kb_nSpace],h_phi);
	      
	      //eps_rho = epsFact_rho*h_phi;
	      //eps_mu  = epsFact_mu *h_phi;

	      eps_rho = epsFact_rho*elementDiameter[eN];
	      eps_mu  = epsFact_mu *elementDiameter[eN];
	      
	      //compute shape and solution information
	      //shape
	      ck.gradTrialFromRef(&  p_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,p_grad_trial_trace);
	      ck.gradTrialFromRef(&vel_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,vel_grad_trial_trace);
	      //ck.gradTrialFromRef(&vel_grad_test_trace_ref [ebN_local_kb_nSpace*nDOF_test_element] ,jacInv_ext,vel_grad_test_trace);
	      //solution and gradients	
	      ck.valFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],&p_trial_trace_ref[ebN_local_kb*nDOF_test_element],p_ext);
	      ck.valFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext);
	      ck.valFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],v_ext);
	      ck.valFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],w_ext);
	      ck.gradFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],p_grad_trial_trace,grad_p_ext);
	      ck.gradFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_u_ext);
	      ck.gradFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_v_ext);
	      ck.gradFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_w_ext);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  int j_nSpace = j*nSpace; 
		  p_test_dS[j] = p_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
		  vel_test_dS[j] = vel_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
		  
		  vel_grad_test_dS[j_nSpace+0] = vel_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_test_element+j_nSpace+0]*dS;  // ido assumes bubnov galerkin
		  vel_grad_test_dS[j_nSpace+1] = vel_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_test_element+j_nSpace+1]*dS;  // ido assumes bubnov galerkin
		  vel_grad_test_dS[j_nSpace+2] = vel_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_test_element+j_nSpace+2]*dS;  // ido assumes bubnov galerkin

		}

	      //
	      //load the boundary values
	      //
	      //isDOFBoundary_p[ebNE_kb] ebqe_bc_p_ext[ebNE_kb];
	      //isDOFBoundary_u[ebNE_kb] ebqe_bc_u_ext[ebNE_kb];
	      //isDOFBoundary_v[ebNE_kb] ebqe_bc_v_ext[ebNE_kb];
	      //isDOFBoundary_w[ebNE_kb] ebqe_bc_w_ext[ebNE_kb];
	      // For now: hard wired zero weak BC 
	      
	      
	      // 
	      //calculate the pde coefficients using the solution and the boundary values for the solution 
	      // 
              gi[0] = 0.0;
	      gi[1] = 0.0;
	      gi[2] = 0.0;	
	
	      unormal = normal[0]*u_ext + normal[1]*v_ext + normal[2]*w_ext;
	      gnormal = normal[0]*gi[0] + normal[1]*gi[1] + normal[2]*gi[2];
	      uneg = 0.5*(unormal - fabs(unormal));
	  
	     	      
	
	      // 
	      //calculate the numerical fluxes 
	      // 

	      //
	      ck.calculateGScale(G,normal,h_penalty);
              penalty =  10.0*eps_mu/h_penalty;
   
	      //
	      //update residuals
	      //
	      for (int i=0;i<nDOF_test_element;i++)
		{
		  int i_nSpace = i*nSpace;
		  elementResidual_p[i] += (unormal - gnormal)*p_test_dS[i];
		  /*
 	          elementResidual_u[i] += (-eps_mu*(grad_u_ext[i_nSpace+0]*normal[0] + grad_u_ext[i_nSpace+0]*normal[0]+ 
		                                    grad_u_ext[i_nSpace+1]*normal[1] + grad_v_ext[i_nSpace+0]*normal[1]+
						    grad_u_ext[i_nSpace+2]*normal[2] + grad_w_ext[i_nSpace+0]*normal[2]) 
						  +penalty*(u_ext-gi[0]) + p_ext*normal[0] - uneg*eps_rho*(u_ext-gi[0]))*vel_test_dS[i]
					+ eps_mu*(  ((v_ext - gi[1])*normal[0] - (u_ext - gi[0])*normal[1])*vel_grad_test_dS[i_nSpace+1] +
					            ((w_ext - gi[2])*normal[0] - (u_ext - gi[0])*normal[2])*vel_grad_test_dS[i_nSpace+2]);	  
						    						  
		  elementResidual_v[i] += (-eps_mu*(grad_v_ext[i_nSpace+0]*normal[0] + grad_u_ext[i_nSpace+1]*normal[0]+ 
		                                    grad_v_ext[i_nSpace+1]*normal[1] + grad_v_ext[i_nSpace+1]*normal[1]+
						    grad_v_ext[i_nSpace+2]*normal[2] + grad_w_ext[i_nSpace+1]*normal[2]) 
						  +penalty*(v_ext-gi[1]) + p_ext*normal[1] - uneg*eps_rho*(v_ext-gi[1]))*vel_test_dS[i]
					+ eps_mu*(  ((u_ext - gi[0])*normal[1] - (v_ext - gi[1])*normal[0])*vel_grad_test_dS[i_nSpace+0] +
					            ((w_ext - gi[2])*normal[1] - (v_ext - gi[1])*normal[2])*vel_grad_test_dS[i_nSpace+2]);
	       
		  elementResidual_w[i] +=(-eps_mu*(grad_w_ext[i_nSpace+0]*normal[0] + grad_u_ext[i_nSpace+2]*normal[0]+ 
		                                   grad_w_ext[i_nSpace+1]*normal[1] + grad_v_ext[i_nSpace+2]*normal[1]+
						   grad_w_ext[i_nSpace+2]*normal[2] + grad_w_ext[i_nSpace+2]*normal[2]) 
						  +penalty*(w_ext-gi[2]) + p_ext*normal[2] - uneg*eps_rho*(w_ext-gi[2]))*vel_test_dS[i]
					+ eps_mu*(  ((u_ext - gi[0])*normal[2] - (w_ext - gi[2])*normal[0])*vel_grad_test_dS[i_nSpace+0] +
					            ((v_ext - gi[1])*normal[2] - (w_ext - gi[2])*normal[1])*vel_grad_test_dS[i_nSpace+1]);*/
		}//i
		
		//penalty = 1.0e-1;
		
		for (int i=0;i<nDOF_test_element;i++)
		{		  
 	          //elementResidual_u[i] += penalty*(u_ext-gi[0])*vel_test_dS[i];  						    						  
		  //elementResidual_v[i] += penalty*(v_ext-gi[1])*vel_test_dS[i];  	       
		  //elementResidual_w[i] += penalty*(w_ext-gi[2])*vel_test_dS[i];  
		}
		
		
	    }//kb
	  //
	  //update the element and global residual storage
	  //
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      int eN_i = eN*nDOF_test_element+i;
	      
	      globalResidual[offset_p+stride_p*p_l2g[eN_i]]  +=elementResidual_p[i];
	      globalResidual[offset_u+stride_u*vel_l2g[eN_i]]+=elementResidual_u[i];
	      globalResidual[offset_v+stride_v*vel_l2g[eN_i]]+=elementResidual_v[i];
	      globalResidual[offset_w+stride_w*vel_l2g[eN_i]]+=elementResidual_w[i];
	    } //i

	}//ebNE
	
    }

    void calculateJacobian(//element
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
			   //physics
			   double* elementDiameter,
			   double hFactor,
			   int nElements_global,
			   double useRBLES,
			   double useMetrics, 
			   double alphaBDF,
			   double epsFact_rho,
			   double epsFact_mu,
			   double sigma,
			   double rho_0,
			   double nu_0,
			   double rho_1,
			   double nu_1,
			   double Ct_sge,
			   double Cd_sge,
			   double C_dc,
			   int* p_l2g, 
			   int* vel_l2g,
			   double* p_dof, double* u_dof, double* v_dof, double* w_dof,
			   double* g,
			   double* phi,
			   double* normal_phi,
			   double* kappa_phi,
			   double* q_mom_u_acc_beta_bdf, double* q_mom_v_acc_beta_bdf, double* q_mom_w_acc_beta_bdf,
			   double* q_velocity_sge,
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
			   double* ebqe_normal_phi_ext,
			   double* ebqe_kappa_phi_ext,
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
			   int* csrColumnOffsets_eb_w_w)
    {

  //
  //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
  //
  for(int eN=0;eN<nElements_global;eN++)
    {
      register double eps_rho,eps_mu;

      register double  elementJacobian_p_p[nDOF_test_element][nDOF_trial_element],
	elementJacobian_p_u[nDOF_test_element][nDOF_trial_element],
	elementJacobian_p_v[nDOF_test_element][nDOF_trial_element],
	elementJacobian_p_w[nDOF_test_element][nDOF_trial_element],
	elementJacobian_u_p[nDOF_test_element][nDOF_trial_element],
	elementJacobian_u_u[nDOF_test_element][nDOF_trial_element],
	elementJacobian_u_v[nDOF_test_element][nDOF_trial_element],
	elementJacobian_u_w[nDOF_test_element][nDOF_trial_element],
	elementJacobian_v_p[nDOF_test_element][nDOF_trial_element],
	elementJacobian_v_u[nDOF_test_element][nDOF_trial_element],
	elementJacobian_v_v[nDOF_test_element][nDOF_trial_element],
	elementJacobian_v_w[nDOF_test_element][nDOF_trial_element],
	elementJacobian_w_p[nDOF_test_element][nDOF_trial_element],
	elementJacobian_w_u[nDOF_test_element][nDOF_trial_element],
	elementJacobian_w_v[nDOF_test_element][nDOF_trial_element],
	elementJacobian_w_w[nDOF_test_element][nDOF_trial_element];
      for (int i=0;i<nDOF_test_element;i++)
	for (int j=0;j<nDOF_trial_element;j++)
	  {
	    elementJacobian_p_p[i][j]=0.0;
	    elementJacobian_p_u[i][j]=0.0;
	    elementJacobian_p_v[i][j]=0.0;
	    elementJacobian_p_w[i][j]=0.0;
	    elementJacobian_u_p[i][j]=0.0;
	    elementJacobian_u_u[i][j]=0.0;
	    elementJacobian_u_v[i][j]=0.0;
	    elementJacobian_u_w[i][j]=0.0;
	    elementJacobian_v_p[i][j]=0.0;
	    elementJacobian_v_u[i][j]=0.0;
	    elementJacobian_v_v[i][j]=0.0;
	    elementJacobian_v_w[i][j]=0.0;
	    elementJacobian_w_p[i][j]=0.0;
	    elementJacobian_w_u[i][j]=0.0;
	    elementJacobian_w_v[i][j]=0.0;
	    elementJacobian_w_w[i][j]=0.0;
	  }
      for  (int k=0;k<nQuadraturePoints_element;k++)
        {
	  int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
	    eN_k_nSpace = eN_k*nSpace,
	    eN_nDOF_trial_element = eN*nDOF_trial_element; //index to a vector at a quadrature point

	  //declare local storage
	  register double p=0.0,u=0.0,v=0.0,w=0.0,
	    grad_p[nSpace],grad_u[nSpace],grad_v[nSpace],grad_w[nSpace],
	    mom_u_acc=0.0,
	    dmom_u_acc_u=0.0,
	    mom_v_acc=0.0,
	    dmom_v_acc_v=0.0,
	    mom_w_acc=0.0,
	    dmom_w_acc_w=0.0,
	    mom_u_acc_t=0.0,
	    dmom_u_acc_u_t=0.0,
	    mom_v_acc_t=0.0,
	    dmom_v_acc_v_t=0.0,
	    mom_w_acc_t=0.0,
	    dmom_w_acc_w_t=0.0,
	    pdeResidual_p=0.0,
	    pdeResidual_u=0.0,
	    pdeResidual_v=0.0,
	    pdeResidual_w=0.0,	    
	    dpdeResidual_p_u[nDOF_trial_element],dpdeResidual_p_v[nDOF_trial_element],dpdeResidual_p_w[nDOF_trial_element],
	    dpdeResidual_u_p[nDOF_trial_element],dpdeResidual_u_u[nDOF_trial_element],
	    dpdeResidual_v_p[nDOF_trial_element],dpdeResidual_v_v[nDOF_trial_element],
	    dpdeResidual_w_p[nDOF_trial_element],dpdeResidual_w_w[nDOF_trial_element],
	    subgrid_p=0.0,
	    subgrid_u=0.0,
	    subgrid_v=0.0,
	    subgrid_w=0.0,	    
	    dsubgrid_p_u[nDOF_trial_element],
	    dsubgrid_p_v[nDOF_trial_element],
	    dsubgrid_p_w[nDOF_trial_element],
	    dsubgrid_u_p[nDOF_trial_element],
	    dsubgrid_u_u[nDOF_trial_element],
	    dsubgrid_v_p[nDOF_trial_element],
	    dsubgrid_v_v[nDOF_trial_element],
	    dsubgrid_w_p[nDOF_trial_element],
	    dsubgrid_w_w[nDOF_trial_element],
	    tau_0=0.0,
	    tau_1=0.0,
	    jac[nSpace*nSpace],
	    jacDet,
	    jacInv[nSpace*nSpace],
	    p_grad_trial[nDOF_trial_element*nSpace],vel_grad_trial[nDOF_trial_element*nSpace],
	    dV,
	    p_test_dV[nDOF_test_element],vel_test_dV[nDOF_test_element],
	    p_grad_test_dV[nDOF_test_element*nSpace],vel_grad_test_dV[nDOF_test_element*nSpace],
	    x,y,z,
	    G[nSpace*nSpace],G_dd_G,tr_G,h_phi,
            rho,mu,H_rho,d_rho,H_mu,d_mu,norm_Rv;
	  //get jacobian, etc for mapping reference element
	  ck.calculateMapping_element(eN,
				      k,
				      mesh_dof,
				      mesh_l2g,
				      mesh_trial_ref,
				      mesh_grad_trial_ref,
				      jac,
				      jacDet,
				      jacInv,
				      x,y,z);
	  //get the physical integration weight
	  dV = fabs(jacDet)*dV_ref[k];
	  ck.calculateG(jacInv,G,G_dd_G,tr_G);
	  ck.calculateGScale(G,&normal_phi[eN_k_nSpace],h_phi);
		 
	  eps_rho = epsFact_rho*h_phi;
	  eps_mu = epsFact_mu*h_phi;

          eps_rho = epsFact_rho*elementDiameter[eN];
	  eps_mu  = epsFact_mu *elementDiameter[eN];
	  
	  //get the trial function gradients
	  ck.gradTrialFromRef(&p_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,p_grad_trial);
	  ck.gradTrialFromRef(&vel_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,vel_grad_trial);
	  //get the solution 	
	  ck.valFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],&p_trial_ref[k*nDOF_trial_element],p);
	  ck.valFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],u);
	  ck.valFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],v);
	  ck.valFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],w);
	  //get the solution gradients
	  ck.gradFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],p_grad_trial,grad_p);
	  ck.gradFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_u);
	  ck.gradFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_v);
	  ck.gradFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_w);
	  //precalculate test function products with integration weights
	  for (int j=0;j<nDOF_trial_element;j++)
	    {
	      p_test_dV[j] = p_test_ref[k*nDOF_trial_element+j]*dV;
	      vel_test_dV[j] = vel_test_ref[k*nDOF_trial_element+j]*dV;
	      for (int I=0;I<nSpace;I++)
		{
		  p_grad_test_dV[j*nSpace+I]   = p_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
		  vel_grad_test_dV[j*nSpace+I] = vel_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin}
		}
	    }


          H_rho = RBLES::smoothedHeaviside(eps_rho,phi[eN_k]);
          d_rho = RBLES::smoothedDirac(eps_rho,phi[eN_k]);
          H_mu = RBLES::smoothedHeaviside(eps_mu,phi[eN_k]);
          d_mu = RBLES::smoothedDirac(eps_mu,phi[eN_k]);
  
          rho = rho_0*(1.0-H_rho)+rho_1*H_rho;
          mu  = nu_0*(1.0-H_mu)+nu_1*H_mu;
	  //cek debug

           //u momentum accumulation
          mom_u_acc=u;
          dmom_u_acc_u=1.0;
  
          //v momentum accumulation
          mom_v_acc=v;
          dmom_v_acc_v=1.0;
  
          //w momentum accumulation
          mom_w_acc=w;
          dmom_w_acc_w=1.0;

          //
          //calculate time derivatives
          //
          ck.bdf(alphaBDF,
		 q_mom_u_acc_beta_bdf[eN_k],
		 mom_u_acc,
		 dmom_u_acc_u,
		 mom_u_acc_t,
		 dmom_u_acc_u_t);
          ck.bdf(alphaBDF,
		 q_mom_v_acc_beta_bdf[eN_k],
		 mom_v_acc,
		 dmom_v_acc_v,
		 mom_v_acc_t,
		 dmom_v_acc_v_t);
          ck.bdf(alphaBDF,
		 q_mom_w_acc_beta_bdf[eN_k],
		 mom_w_acc,
		 dmom_w_acc_w,
		 mom_w_acc_t,
		 dmom_w_acc_w_t);
          //
          //calculate subgrid error contribution to the Jacobian (strong residual, adjoint, jacobian of strong residual)
          //
//u = v = w = 1.0;
/*
		 mom_u_acc= 0.0;
		 dmom_u_acc_u= 0.0;
		 mom_u_acc_t= 0.0;
		 dmom_u_acc_u_t= 0.0;
		 mom_v_acc= 0.0;
		 dmom_v_acc_v= 0.0;
		 mom_v_acc_t= 0.0;
		 dmom_v_acc_v_t= 0.0;
        
		 mom_w_acc= 0.0;
		 dmom_w_acc_w= 0.0;
		 mom_w_acc_t= 0.0;
		 dmom_w_acc_w_t= 0.0;*/



          //calculate tau and tau*Res					      
          double Dt=(dmom_u_acc_u_t/dmom_u_acc_u);
          double v_d_Gv= 0.0;
	  double vel[nSpace] = {u,v,w};
          for(int I=0;I<nSpace;I++)
            for (int J=0;J<nSpace;J++)
	      v_d_Gv += vel[I]*G[I*nSpace+J]*vel[J];

    
          tau_0 = 1.0/sqrt(Ct_sge*rho*rho*Dt*Dt + rho*rho*v_d_Gv + Cd_sge*mu*mu*G_dd_G);
          tau_0 = 1.0/sqrt(Ct_sge*rho*rho*Dt*Dt + Cd_sge*mu*mu*G_dd_G);
          tau_1 = 1.0/(tr_G*tau_0);
	  //tau_0 = 0.0;      
          //tau_1 = 0.0; 		
		
          //
          //calculate strong residual
          //
	  pdeResidual_p = grad_u[0] + grad_v[1] + grad_w[2];

	  pdeResidual_u = rho *(mom_u_acc_t + u*grad_u[0] + v*grad_u[1] + w*grad_u[2] - g[0])
	                + grad_p[0]; // - mu poisson_u  ==> ommited as only first order right now !! 

	  pdeResidual_v = rho *(mom_v_acc_t + u*grad_v[0] + v*grad_v[1] + w*grad_v[2] - g[1])
	                + grad_p[1]; // - mu poisson_u  ==> ommited as only first order right now !! 
	                
	  pdeResidual_w = rho *(mom_w_acc_t + u*grad_w[0] + v*grad_w[1] + w*grad_w[2] - g[2])
	                + grad_p[2]; // - mu poisson_u  ==> ommited as only first order right now !! 


	  pdeResidual_u = rho *(mom_u_acc_t - g[0])
	                + grad_p[0]; // - mu poisson_u  ==> ommited as only first order right now !! 

	  pdeResidual_v = rho *(mom_v_acc_t - g[1])
	                + grad_p[1]; // - mu poisson_u  ==> ommited as only first order right now !! 
	                
	  pdeResidual_w = rho *(mom_w_acc_t - g[2])
	                + grad_p[2]; // - mu poisson_u  ==> ommited as only first order right now !! 

          //
          //calculate reconstructed subgrid scales
          //	                        
          subgrid_u = -tau_0*pdeResidual_u;
          subgrid_v = -tau_0*pdeResidual_v;
          subgrid_w = -tau_0*pdeResidual_w;

          subgrid_p = -tau_1*pdeResidual_p;
 
          
          //calculate Jacobian of reconstructed subgrid scales
          for (int j=0;j<nDOF_trial_element;j++)
            {
	      int j_nSpace = j*nSpace;

	      dpdeResidual_p_u[j] = vel_grad_trial[j_nSpace+0];
	      dpdeResidual_p_v[j] = vel_grad_trial[j_nSpace+1];
	      dpdeResidual_p_w[j] = vel_grad_trial[j_nSpace+2];

	      dpdeResidual_u_p[j]=p_grad_trial[j_nSpace+0];
	      dpdeResidual_u_u[j]=rho*(dmom_u_acc_u_t*vel_trial_ref[j]
	                            );// + u*vel_grad_trial[j_nSpace+0] 
	                             //+ v*vel_grad_trial[j_nSpace+1] 
	                             //+ w*vel_grad_trial[j_nSpace+2]);

	      dpdeResidual_v_p[j]=p_grad_trial[j_nSpace+1];
	      dpdeResidual_v_v[j]=rho*(dmom_v_acc_v_t*vel_trial_ref[j]
	                          );//   + u*vel_grad_trial[j_nSpace+0] 
	                          //   + v*vel_grad_trial[j_nSpace+1] 
	                          //   + w*vel_grad_trial[j_nSpace+2]);
	                             
	      dpdeResidual_w_p[j]=p_grad_trial[j_nSpace+2];
	      dpdeResidual_w_w[j]=rho*(dmom_w_acc_w_t*vel_trial_ref[j]
	                           );//  + u*vel_grad_trial[j_nSpace+0] 
	                           //  + v*vel_grad_trial[j_nSpace+1] 
	                           //  + w*vel_grad_trial[j_nSpace+2]);
              /* p */
              dsubgrid_p_u[j] = -tau_1*dpdeResidual_p_u[j];
              dsubgrid_p_v[j] = -tau_1*dpdeResidual_p_v[j];
              dsubgrid_p_w[j] = -tau_1*dpdeResidual_p_w[j];

              /* u */
              dsubgrid_u_p[j] = -tau_0*dpdeResidual_u_p[j];
              dsubgrid_u_u[j] = -tau_0*dpdeResidual_u_u[j];
              /* v */
              dsubgrid_v_p[j] = -tau_0*dpdeResidual_v_p[j];
              dsubgrid_v_v[j] = -tau_0*dpdeResidual_v_v[j];
              /* w */
              dsubgrid_w_p[j] = -tau_0*dpdeResidual_w_p[j];
              dsubgrid_w_w[j] = -tau_0*dpdeResidual_w_w[j];


           }
           //
          //calculate discontinuity capturing addition to viscosity
          //	   
	  norm_Rv = sqrt(pdeResidual_u*pdeResidual_u + pdeResidual_v*pdeResidual_v + pdeResidual_w*pdeResidual_w);
//	  mu += C_dc*norm_Rv/sqrt(G_dd_G);


	  //cek todo add RBLES terms consistent to residual modifications or ignore them partials w.r.t the additional RBLES terms
  	  for(int i=0;i<nDOF_test_element;i++)
	    {
	      int i_nSpace=i*nSpace;
	      for(int j=0;j<nDOF_trial_element;j++) 
		{ 
		  int j_nSpace = j*nSpace;
/*
		  elementJacobian_u_u[i][j] += vel_trial_ref[j]*vel_test_dV[i];
		  elementJacobian_v_v[i][j] += vel_trial_ref[j]*vel_test_dV[i];
		  elementJacobian_w_w[i][j] += vel_trial_ref[j]*vel_test_dV[i];


		  elementJacobian_p_p[i][j] += p_grad_trial[j_nSpace+0]*p_grad_test_dV[i_nSpace+0]
	                                      +p_grad_trial[j_nSpace+1]*p_grad_test_dV[i_nSpace+1]
	                                      +p_grad_trial[j_nSpace+2]*p_grad_test_dV[i_nSpace+2] + p_trial_ref[j]*p_test_dV[i];
*/					       

          // p 
		  elementJacobian_p_p[i][j] += dsubgrid_u_p[j]*p_grad_test_dV[i_nSpace+0]
	                                      +dsubgrid_v_p[j]*p_grad_test_dV[i_nSpace+1]
	                                      +dsubgrid_w_p[j]*p_grad_test_dV[i_nSpace+2];
					       
     
		  elementJacobian_p_u[i][j] += - vel_grad_trial[j_nSpace+0]*p_test_dV[i]
		                               + dsubgrid_u_u[j]*p_grad_test_dV[i_nSpace+0];
		  elementJacobian_p_v[i][j] += - vel_grad_trial[j_nSpace+1]*p_test_dV[i]
		                               + dsubgrid_v_v[j]*p_grad_test_dV[i_nSpace+1];
		  elementJacobian_p_w[i][j] += - vel_grad_trial[j_nSpace+2]*p_test_dV[i]
		                               + dsubgrid_w_w[j]*p_grad_test_dV[i_nSpace+2];



         // u 
		  elementJacobian_u_u[i][j] += dpdeResidual_u_u[j]*vel_test_dV[i]	 
		          ;//              -rho*vel_grad_test_dV[i_nSpace+0] * (dsubgrid_u_u[j]*u)// + vel_trial_ref[j]*subgrid_u)
		            //            -rho*vel_grad_test_dV[i_nSpace+1] * (dsubgrid_u_u[j]*v)// + vel_trial_ref[j]*subgrid_v)
		             //           -rho*vel_grad_test_dV[i_nSpace+2] * (dsubgrid_u_u[j]*w);// + vel_trial_ref[j]*subgrid_w);

		  elementJacobian_u_p[i][j] += -p_trial_ref[j]*vel_grad_test_dV[i_nSpace+0]	                                  
		               ;//         -rho*vel_grad_test_dV[i_nSpace+0] * (dsubgrid_u_p[j]*u)// + v*dsubgrid_u_p[j])
		                //        -rho*vel_grad_test_dV[i_nSpace+1] * (dsubgrid_u_p[j]*v)//+ v*dsubgrid_v_p[j])
		                //        -rho*vel_grad_test_dV[i_nSpace+2] * (dsubgrid_u_p[j]*w);// + v*dsubgrid_w_p[j]);

		  elementJacobian_u_u[i][j] += -dsubgrid_p_u[j]*vel_grad_test_dV[i_nSpace+0]; 		  
		  elementJacobian_u_v[i][j] += -dsubgrid_p_v[j]*vel_grad_test_dV[i_nSpace+0];
		  elementJacobian_u_w[i][j] += -dsubgrid_p_w[j]*vel_grad_test_dV[i_nSpace+0];		

		  elementJacobian_u_u[i][j] += 2.0*mu*vel_grad_trial[j_nSpace+0]*vel_grad_test_dV[i_nSpace+0]
		                                +  mu*vel_grad_trial[j_nSpace+1]*vel_grad_test_dV[i_nSpace+1]			   
		                                +  mu*vel_grad_trial[j_nSpace+2]*vel_grad_test_dV[i_nSpace+2];

		  elementJacobian_u_v[i][j] += mu*vel_grad_trial[j_nSpace+0]*vel_grad_test_dV[i_nSpace+1];				   
		  elementJacobian_u_w[i][j] += mu*vel_grad_trial[j_nSpace+0]*vel_grad_test_dV[i_nSpace+2];


		  


          // v 
		  elementJacobian_v_v[i][j] += dpdeResidual_v_v[j]*vel_test_dV[i]	 
		                   ;//     -rho*vel_grad_test_dV[i_nSpace+0] * (dsubgrid_v_v[j]*u)// + vel_trial_ref[j]*subgrid_u)
		                    //    -rho*vel_grad_test_dV[i_nSpace+1] * (dsubgrid_v_v[j]*v)// + vel_trial_ref[j]*subgrid_v)
		                     //   -rho*vel_grad_test_dV[i_nSpace+2] * (dsubgrid_v_v[j]*w);// + vel_trial_ref[j]*subgrid_w);
					
		  elementJacobian_v_p[i][j] += -p_trial_ref[j]*vel_grad_test_dV[i_nSpace+1]	                                  
		                   ;//     -rho*vel_grad_test_dV[i_nSpace+0] * (dsubgrid_v_p[j]*u)// + v*dsubgrid_u_p[j])
		                   //     -rho*vel_grad_test_dV[i_nSpace+1] * (dsubgrid_v_p[j]*v)// + v*dsubgrid_v_p[j])
		                   //     -rho*vel_grad_test_dV[i_nSpace+2] * (dsubgrid_v_p[j]*w);// + v*dsubgrid_w_p[j]);

		  elementJacobian_v_u[i][j] += -dsubgrid_p_u[j]*vel_grad_test_dV[i_nSpace+1]; 		  
		  elementJacobian_v_v[i][j] += -dsubgrid_p_v[j]*vel_grad_test_dV[i_nSpace+1];
		  elementJacobian_v_w[i][j] += -dsubgrid_p_w[j]*vel_grad_test_dV[i_nSpace+1];		  

		  elementJacobian_v_u[i][j] += mu*vel_grad_trial[j_nSpace+1]*vel_grad_test_dV[i_nSpace+0];
		  elementJacobian_v_v[i][j] +=     mu*vel_grad_trial[j_nSpace+0]*vel_grad_test_dV[i_nSpace+0]
		                             + 2.0*mu*vel_grad_trial[j_nSpace+1]*vel_grad_test_dV[i_nSpace+1]				   
		                                +  mu*vel_grad_trial[j_nSpace+2]*vel_grad_test_dV[i_nSpace+2];				   
		  elementJacobian_v_w[i][j] += mu*vel_grad_trial[j_nSpace+1]*vel_grad_test_dV[i_nSpace+2];
		  		  


		 
          // w 
		  elementJacobian_w_w[i][j] += dpdeResidual_w_w[j]*vel_test_dV[i]		 
		                     ;//   -rho*vel_grad_test_dV[i_nSpace+0] * (dsubgrid_w_w[j]*u)// + vel_trial_ref[j]*subgrid_u)
		                      //  -rho*vel_grad_test_dV[i_nSpace+1] * (dsubgrid_w_w[j]*v)// + vel_trial_ref[j]*subgrid_v)
		                     //   -rho*vel_grad_test_dV[i_nSpace+2] * (dsubgrid_w_w[j]*w);// + vel_trial_ref[j]*subgrid_w);
		  
		  elementJacobian_w_p[i][j] += -p_trial_ref[j]*vel_grad_test_dV[i_nSpace+2]	                                  
		                      ;//  -rho*vel_grad_test_dV[i_nSpace+0] * (dsubgrid_w_p[j]*u)// + v*dsubgrid_u_p[j])
		                       // -rho*vel_grad_test_dV[i_nSpace+1] * (dsubgrid_w_p[j]*v)// + v*dsubgrid_v_p[j])
		                       // -rho*vel_grad_test_dV[i_nSpace+2] * (dsubgrid_w_p[j]*w);// + v*dsubgrid_w_p[j]);

		  elementJacobian_w_u[i][j] += -dsubgrid_p_u[j]*vel_grad_test_dV[i_nSpace+2];
		  elementJacobian_w_v[i][j] += -dsubgrid_p_v[j]*vel_grad_test_dV[i_nSpace+2];
		  elementJacobian_w_w[i][j] += -dsubgrid_p_w[j]*vel_grad_test_dV[i_nSpace+2];

		  elementJacobian_w_u[i][j] += mu*vel_grad_trial[j_nSpace+2]*vel_grad_test_dV[i_nSpace+0];				   
		  elementJacobian_w_v[i][j] += mu*vel_grad_trial[j_nSpace+2]*vel_grad_test_dV[i_nSpace+1];
		  
		  elementJacobian_w_w[i][j] +=     mu*vel_grad_trial[j_nSpace+0]*vel_grad_test_dV[i_nSpace+0]
		                                +  mu*vel_grad_trial[j_nSpace+1]*vel_grad_test_dV[i_nSpace+1]				   
		                             + 2.0*mu*vel_grad_trial[j_nSpace+2]*vel_grad_test_dV[i_nSpace+2];


		


		}//j
            }//i
	}//k
	  //
	  //load into element Jacobian into global Jacobian
	  //
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      register int eN_i = eN*nDOF_test_element+i;
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  register int eN_i_j = eN_i*nDOF_trial_element+j;
		  globalJacobian[csrRowIndeces_p_p[eN_i] + csrColumnOffsets_p_p[eN_i_j]] += elementJacobian_p_p[i][j];
		  globalJacobian[csrRowIndeces_p_u[eN_i] + csrColumnOffsets_p_u[eN_i_j]] += elementJacobian_p_u[i][j];
		  globalJacobian[csrRowIndeces_p_v[eN_i] + csrColumnOffsets_p_v[eN_i_j]] += elementJacobian_p_v[i][j];
		  globalJacobian[csrRowIndeces_p_w[eN_i] + csrColumnOffsets_p_w[eN_i_j]] += elementJacobian_p_w[i][j];

		  globalJacobian[csrRowIndeces_u_p[eN_i] + csrColumnOffsets_u_p[eN_i_j]] += elementJacobian_u_p[i][j];
		  globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] += elementJacobian_u_u[i][j];
		  globalJacobian[csrRowIndeces_u_v[eN_i] + csrColumnOffsets_u_v[eN_i_j]] += elementJacobian_u_v[i][j];
		  globalJacobian[csrRowIndeces_u_w[eN_i] + csrColumnOffsets_u_w[eN_i_j]] += elementJacobian_u_w[i][j];

		  globalJacobian[csrRowIndeces_v_p[eN_i] + csrColumnOffsets_v_p[eN_i_j]] += elementJacobian_v_p[i][j];
		  globalJacobian[csrRowIndeces_v_u[eN_i] + csrColumnOffsets_v_u[eN_i_j]] += elementJacobian_v_u[i][j];
		  globalJacobian[csrRowIndeces_v_v[eN_i] + csrColumnOffsets_v_v[eN_i_j]] += elementJacobian_v_v[i][j];
		  globalJacobian[csrRowIndeces_v_w[eN_i] + csrColumnOffsets_v_w[eN_i_j]] += elementJacobian_v_w[i][j];

		  globalJacobian[csrRowIndeces_w_p[eN_i] + csrColumnOffsets_w_p[eN_i_j]] += elementJacobian_w_p[i][j];
		  globalJacobian[csrRowIndeces_w_u[eN_i] + csrColumnOffsets_w_u[eN_i_j]] += elementJacobian_w_u[i][j];
		  globalJacobian[csrRowIndeces_w_v[eN_i] + csrColumnOffsets_w_v[eN_i_j]] += elementJacobian_w_v[i][j];
		  globalJacobian[csrRowIndeces_w_w[eN_i] + csrColumnOffsets_w_w[eN_i_j]] += elementJacobian_w_w[i][j];
		}//j
	    }//i
	}//elements
      //
      //loop over exterior element boundaries to compute the surface integrals and load them into the global Jacobian
      //
      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
	{ 
	  register int ebN = exteriorElementBoundariesArray[ebNE],
	    eN  = elementBoundaryElementsArray[ebN*2+0],
	    eN_nDOF_trial_element = eN*nDOF_trial_element,
	    ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0];
	  register double eps_rho,eps_mu;

	 double elementJacobian_p_p[nDOF_test_element][nDOF_trial_element],
		elementJacobian_p_u[nDOF_test_element][nDOF_trial_element],
		elementJacobian_p_v[nDOF_test_element][nDOF_trial_element],
		elementJacobian_p_w[nDOF_test_element][nDOF_trial_element],
		elementJacobian_u_p[nDOF_test_element][nDOF_trial_element],
		elementJacobian_u_u[nDOF_test_element][nDOF_trial_element],
		elementJacobian_u_v[nDOF_test_element][nDOF_trial_element],
		elementJacobian_u_w[nDOF_test_element][nDOF_trial_element],
		elementJacobian_v_p[nDOF_test_element][nDOF_trial_element],
		elementJacobian_v_u[nDOF_test_element][nDOF_trial_element],
		elementJacobian_v_v[nDOF_test_element][nDOF_trial_element],
		elementJacobian_v_w[nDOF_test_element][nDOF_trial_element],
		elementJacobian_w_p[nDOF_test_element][nDOF_trial_element],
		elementJacobian_w_u[nDOF_test_element][nDOF_trial_element],
		elementJacobian_w_v[nDOF_test_element][nDOF_trial_element],
		elementJacobian_w_w[nDOF_test_element][nDOF_trial_element];
	
      for (int i=0;i<nDOF_test_element;i++)
	for (int j=0;j<nDOF_trial_element;j++)
	  {
	    elementJacobian_p_p[i][j]=0.0;
	    elementJacobian_p_u[i][j]=0.0;
	    elementJacobian_p_v[i][j]=0.0;
	    elementJacobian_p_w[i][j]=0.0;
	    elementJacobian_u_p[i][j]=0.0;
	    elementJacobian_u_u[i][j]=0.0;
	    elementJacobian_u_v[i][j]=0.0;
	    elementJacobian_u_w[i][j]=0.0;
	    elementJacobian_v_p[i][j]=0.0;
	    elementJacobian_v_u[i][j]=0.0;
	    elementJacobian_v_v[i][j]=0.0;
	    elementJacobian_v_w[i][j]=0.0;
	    elementJacobian_w_p[i][j]=0.0;
	    elementJacobian_w_u[i][j]=0.0;
	    elementJacobian_w_v[i][j]=0.0;
	    elementJacobian_w_w[i][j]=0.0;
	  }
	
	
	  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	    { 
	      register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
		ebNE_kb_nSpace = ebNE_kb*nSpace,
		ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
		ebN_local_kb_nSpace = ebN_local_kb*nSpace;

	      register double p_ext=0.0,
		u_ext=0.0,
		v_ext=0.0,
		w_ext=0.0,
		grad_p_ext[nSpace],
		grad_u_ext[nSpace],
		grad_v_ext[nSpace],
		grad_w_ext[nSpace],
		jac_ext[nSpace*nSpace],
		jacDet_ext,
		jacInv_ext[nSpace*nSpace],
		boundaryJac[nSpace*(nSpace-1)],
		metricTensor[(nSpace-1)*(nSpace-1)],
		metricTensorDetSqrt,
		p_grad_trial_trace[nDOF_trial_element*nSpace],
		vel_grad_trial_trace[nDOF_trial_element*nSpace],
		dS,
		p_test_dS[nDOF_test_element],
		vel_test_dS[nDOF_test_element],vel_grad_test_dS[nDOF_test_element*nSpace],
		normal[3],
		x_ext,y_ext,z_ext,
		G[nSpace*nSpace],G_dd_G,tr_G,h_phi,h_penalty,penalty,gi[nSpace], unormal,gnormal,uneg;
	      ck.calculateMapping_elementBoundary(eN,
						  ebN_local,
						  kb,
						  ebN_local_kb,
						  mesh_dof,
						  mesh_l2g,
						  mesh_trial_trace_ref,
						  mesh_grad_trial_trace_ref,
						  boundaryJac_ref,
						  jac_ext,
						  jacDet_ext,
						  jacInv_ext,
						  boundaryJac,
						  metricTensor,
						  metricTensorDetSqrt,
						  normal_ref,
						  normal,
						  x_ext,y_ext,z_ext);
	      dS = metricTensorDetSqrt*dS_ref[kb];
	      ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
	      ck.calculateGScale(G,&ebqe_normal_phi_ext[ebNE_kb_nSpace],h_phi);

	      eps_rho = epsFact_rho*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN]);
	      eps_mu  = epsFact_mu *(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN]);

	      eps_rho = epsFact_rho*elementDiameter[eN];
	      eps_mu  = epsFact_mu *elementDiameter[eN];
	      
	      //compute shape and solution information
	      //shape
	      ck.gradTrialFromRef(&p_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,p_grad_trial_trace);
	      ck.gradTrialFromRef(&vel_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,vel_grad_trial_trace);
	      //solution and gradients	
	      ck.valFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],&p_trial_trace_ref[ebN_local_kb*nDOF_test_element],p_ext);
	      ck.valFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext);
	      ck.valFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],v_ext);
	      ck.valFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],w_ext);
	      ck.gradFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],p_grad_trial_trace,grad_p_ext);
	      ck.gradFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_u_ext);
	      ck.gradFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_v_ext);
	      ck.gradFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_w_ext);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  int j_nSpace = j*nSpace; 
		  p_test_dS[j] = p_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
		  vel_test_dS[j] = vel_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
		  
		  vel_grad_test_dS[j_nSpace+0] = vel_grad_test_trace_ref[ebN_local_kb_nSpace*nDOF_test_element+j_nSpace+0]*dS;  
		  vel_grad_test_dS[j_nSpace+1] = vel_grad_test_trace_ref[ebN_local_kb_nSpace*nDOF_test_element+j_nSpace+1]*dS; 
		  vel_grad_test_dS[j_nSpace+2] = vel_grad_test_trace_ref[ebN_local_kb_nSpace*nDOF_test_element+j_nSpace+2]*dS; 

		}

	      //
	      //load the boundary values
	      //
	      //isDOFBoundary_p[ebNE_kb] ebqe_bc_p_ext[ebNE_kb];
	      //isDOFBoundary_u[ebNE_kb] ebqe_bc_u_ext[ebNE_kb];
	      //isDOFBoundary_v[ebNE_kb] ebqe_bc_v_ext[ebNE_kb];
	      //isDOFBoundary_w[ebNE_kb] ebqe_bc_w_ext[ebNE_kb];
	      // For now: hard wired zero weak BC 
              gi[0] = 0.0;
	      gi[1] = 0.0;
	      gi[2] = 0.0;	      
	      
	      // 
	      //calculate the pde coefficients using the solution and the boundary values for the solution 
	      // 	
	
	      unormal = normal[0]*u_ext + normal[1]*v_ext + normal[2]*w_ext;
	      gnormal = normal[0]*gi[0] + normal[1]*gi[1] + normal[2]*gi[2];
	      uneg = 0.5*(unormal - fabs(unormal));
	     	      
	
	      // 
	      //calculate the numerical fluxes 
	      // 

	      //
	      ck.calculateGScale(G,normal,h_penalty);
              penalty =  10.0*eps_mu/h_penalty;

	      // 
	      //calculate the internal and external trace of the pde coefficients 
	      // 

	     /* for (int i=0;i<nDOF_test_element;i++)
	      for (int i=0;i<nDOF_test_element;i++)
		{
		  int i_nSpace = i*nSpace;
		  elementResidual_p[i] += (unormal - gnormal)*p_test_dS[i];
		  
 	          elementResidual_u[i] += (-eps_mu*(grad_u_ext[i_nSpace+0]*normal[0] + grad_u_ext[i_nSpace+0]*normal[0]+ 
		                                    grad_u_ext[i_nSpace+1]*normal[1] + grad_v_ext[i_nSpace+0]*normal[1]+
						    grad_u_ext[i_nSpace+2]*normal[2] + grad_w_ext[i_nSpace+0]*normal[2]) 
						  +penalty*(u_ext-gi[0]) + p_ext*normal[0] - uneg*eps_rho*(u_ext-gi[0]))*vel_test_dS[i]
					+ eps_mu*(  ((v_ext - gi[1])*normal[0] - (u_ext - gi[0])*normal[1])*vel_grad_test_dS[i_nSpace+1] +
					            ((w_ext - gi[2])*normal[0] - (u_ext - gi[0])*normal[2])*vel_grad_test_dS[i_nSpace+2]);	  
						    						  
		  elementResidual_v[i] += (-eps_mu*(grad_v_ext[i_nSpace+0]*normal[0] + grad_u_ext[i_nSpace+1]*normal[0]+ 
		                                    grad_v_ext[i_nSpace+1]*normal[1] + grad_v_ext[i_nSpace+1]*normal[1]+
						    grad_v_ext[i_nSpace+2]*normal[2] + grad_w_ext[i_nSpace+1]*normal[2]) 
						  +penalty*(v_ext-gi[1]) + p_ext*normal[1] - uneg*eps_rho*(v_ext-gi[1]))*vel_test_dS[i]
					+ eps_mu*(  ((u_ext - gi[0])*normal[1] - (v_ext - gi[1])*normal[0])*vel_grad_test_dS[i_nSpace+0] +
					            ((w_ext - gi[2])*normal[1] - (v_ext - gi[1])*normal[2])*vel_grad_test_dS[i_nSpace+2]);
	       
		  elementResidual_w[i] +=(-eps_mu*(grad_w_ext[i_nSpace+0]*normal[0] + grad_u_ext[i_nSpace+2]*normal[0]+ 
		                                   grad_w_ext[i_nSpace+1]*normal[1] + grad_v_ext[i_nSpace+2]*normal[1]+
						   grad_w_ext[i_nSpace+2]*normal[2] + grad_w_ext[i_nSpace+2]*normal[2]) 
						  +penalty*(w_ext-gi[2]) + p_ext*normal[2] - uneg*eps_rho*(w_ext-gi[2]))*vel_test_dS[i]
					+ eps_mu*(  ((u_ext - gi[0])*normal[2] - (w_ext - gi[2])*normal[0])*vel_grad_test_dS[i_nSpace+0] +
					            ((v_ext - gi[1])*normal[2] - (w_ext - gi[2])*normal[1])*vel_grad_test_dS[i_nSpace+1]);
		}//i
		}//i*/

	      //
	      //calculate the flux jacobian
	      //
	      // vel_trial_trace_ref[ebN_local_kb*nDOF_test_element]
	      
  	    for(int i=0;i<nDOF_test_element;i++)
	      {
	      //int i_nSpace=i*nSpace;
	      for(int j=0;j<nDOF_trial_element;j++) 
		{ 
		  //int j_nSpace = j*nSpace;
		  // //cek debug
		   elementJacobian_p_p[i][j]=0.0;
		   elementJacobian_p_u[i][j]+= normal[0]*vel_trial_trace_ref[ebN_local_kb*nDOF_test_element+j]*p_test_dS[i];
		   elementJacobian_p_v[i][j]+= normal[1]*vel_trial_trace_ref[ebN_local_kb*nDOF_test_element+j]*p_test_dS[i];
		   elementJacobian_p_w[i][j]+= normal[2]*vel_trial_trace_ref[ebN_local_kb*nDOF_test_element+j]*p_test_dS[i];

		  /* elementJacobian_u_p[i][j]+= normal[0]*p_trial_trace_ref[j]*vel_test_dS[i];
		   elementJacobian_u_u[i][j]+= (penalty - uneg*eps_rho)*vel_trial_trace_ref[j]*vel_test_dS[i];
		   elementJacobian_u_v[i][j]=0.0;
		   elementJacobian_u_w[i][j]=0.0;

		   elementJacobian_v_p[i][j]+= normal[1]*p_trial_trace_ref[j]*vel_test_dS[i];
		   elementJacobian_v_u[i][j]=0.0;
		   elementJacobian_v_v[i][j]+= (penalty - uneg*eps_rho)*vel_trial_trace_ref[j]*vel_test_dS[i];
		   elementJacobian_v_w[i][j]=0.0;

		   elementJacobian_w_p[i][j]+= normal[2]*p_trial_trace_ref[j]*vel_test_dS[i];
		   elementJacobian_w_u[i][j]=0.0;
		   elementJacobian_w_v[i][j]=0.0;
		   elementJacobian_w_w[i][j]+= (penalty - uneg*eps_rho)*vel_trial_trace_ref[j]*vel_test_dS[i];*/
		  // //cek debug
		}//j
	      }//i
	      
	      		penalty = 1.0e-1;
		
  	    for(int i=0;i<nDOF_test_element;i++)
	      {
	      //int i_nSpace=i*nSpace;
	      for(int j=0;j<nDOF_trial_element;j++) 
		{ 
		  //int j_nSpace = j*nSpace;
		  // //cek debug


		
		   //elementJacobian_u_u[i][j] += penalty*vel_trial_trace_ref[ebN_local_kb*nDOF_test_element+j]*vel_test_dS[i];
                   //elementJacobian_v_v[i][j] += penalty*vel_trial_trace_ref[ebN_local_kb*nDOF_test_element+j]*vel_test_dS[i];
                   //elementJacobian_w_w[i][j] += penalty*vel_trial_trace_ref[ebN_local_kb*nDOF_test_element+j]*vel_test_dS[i];
		  // //cek debug
		}//j
	      }//i
	      
	      //
	      //update the global Jacobian from the flux Jacobian
	      //
	      
	      
	      
	      for (int i=0;i<nDOF_test_element;i++)
		{
		  register int eN_i = eN*nDOF_test_element+i;
		  for (int j=0;j<nDOF_trial_element;j++)
		    {
		      register int ebN_i_j = ebN*4*nDOF_test_X_trial_element + i*nDOF_trial_element + j;
		  
		      globalJacobian[csrRowIndeces_p_p[eN_i] + csrColumnOffsets_eb_p_p[ebN_i_j]] -=  elementJacobian_p_p[i][j];
		      globalJacobian[csrRowIndeces_p_u[eN_i] + csrColumnOffsets_eb_p_u[ebN_i_j]] -=  elementJacobian_p_u[i][j];
		      globalJacobian[csrRowIndeces_p_v[eN_i] + csrColumnOffsets_eb_p_v[ebN_i_j]] -=  elementJacobian_p_v[i][j];
		      globalJacobian[csrRowIndeces_p_w[eN_i] + csrColumnOffsets_eb_p_w[ebN_i_j]] -=  elementJacobian_p_w[i][j];
		   
		      globalJacobian[csrRowIndeces_u_p[eN_i] + csrColumnOffsets_eb_u_p[ebN_i_j]] +=  elementJacobian_u_p[i][j];
		      globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] +=  elementJacobian_u_u[i][j];
		      globalJacobian[csrRowIndeces_u_v[eN_i] + csrColumnOffsets_eb_u_v[ebN_i_j]] +=  elementJacobian_u_v[i][j];
		      globalJacobian[csrRowIndeces_u_w[eN_i] + csrColumnOffsets_eb_u_w[ebN_i_j]] +=  elementJacobian_u_w[i][j];
		   
		      globalJacobian[csrRowIndeces_v_p[eN_i] + csrColumnOffsets_eb_v_p[ebN_i_j]] +=  elementJacobian_v_p[i][j];
		      globalJacobian[csrRowIndeces_v_u[eN_i] + csrColumnOffsets_eb_v_u[ebN_i_j]] +=  elementJacobian_v_u[i][j];
		      globalJacobian[csrRowIndeces_v_v[eN_i] + csrColumnOffsets_eb_v_v[ebN_i_j]] +=  elementJacobian_v_v[i][j];
		      globalJacobian[csrRowIndeces_v_w[eN_i] + csrColumnOffsets_eb_v_w[ebN_i_j]] +=  elementJacobian_v_w[i][j];
		   
		      globalJacobian[csrRowIndeces_w_p[eN_i] + csrColumnOffsets_eb_w_p[ebN_i_j]] +=  elementJacobian_w_p[i][j];
		      globalJacobian[csrRowIndeces_w_u[eN_i] + csrColumnOffsets_eb_w_u[ebN_i_j]] +=  elementJacobian_w_u[i][j];
		      globalJacobian[csrRowIndeces_w_v[eN_i] + csrColumnOffsets_eb_w_v[ebN_i_j]] +=  elementJacobian_w_v[i][j];
		      globalJacobian[csrRowIndeces_w_w[eN_i] + csrColumnOffsets_eb_w_w[ebN_i_j]] +=  elementJacobian_w_w[i][j];
		    }//j
		}//i
		
		
		
	    }//kb
	}//ebNE
    }//computeJacobian

    void calculateVelocityAverage(int nExteriorElementBoundaries_global,
    				  int* exteriorElementBoundariesArray,
    				  int nInteriorElementBoundaries_global,
    				  int* interiorElementBoundariesArray,
    				  int* elementBoundaryElementsArray,
    				  int* elementBoundaryLocalElementBoundariesArray,
    				  double* mesh_dof,
    				  int* mesh_l2g,
    				  double* mesh_trial_trace_ref,
    				  double* mesh_grad_trial_trace_ref,
    				  double* normal_ref,
    				  double* boundaryJac_ref,
    				  int* vel_l2g,
    				  double* u_dof,
    				  double* v_dof,
    				  double* w_dof,
    				  double* vel_trial_trace_ref,
    				  double* ebqe_velocity,
    				  double* velocityAverage)
    {
      int permutations[nQuadraturePoints_elementBoundary];
      double xArray_left[nQuadraturePoints_elementBoundary*3],
    	xArray_right[nQuadraturePoints_elementBoundary*3];
      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++)
    	{
    	  register int ebN = exteriorElementBoundariesArray[ebNE];
    	  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
    	    {
    	      register int ebN_kb_nSpace = ebN*nQuadraturePoints_elementBoundary*nSpace+kb*nSpace,
    		ebNE_kb_nSpace = ebNE*nQuadraturePoints_elementBoundary*nSpace+kb*nSpace;
    	      velocityAverage[ebN_kb_nSpace+0]=ebqe_velocity[ebNE_kb_nSpace+0];
    	      velocityAverage[ebN_kb_nSpace+1]=ebqe_velocity[ebNE_kb_nSpace+1];
    	      velocityAverage[ebN_kb_nSpace+2]=ebqe_velocity[ebNE_kb_nSpace+2];
    	    }//ebNE
    	}
      for (int ebNI = 0; ebNI < nInteriorElementBoundaries_global; ebNI++)
    	{
    	  register int ebN = interiorElementBoundariesArray[ebNI],
    	    left_eN_global   = elementBoundaryElementsArray[ebN*2+0],
    	    left_ebN_element  = elementBoundaryLocalElementBoundariesArray[ebN*2+0],
    	    right_eN_global  = elementBoundaryElementsArray[ebN*2+1],
    	    right_ebN_element = elementBoundaryLocalElementBoundariesArray[ebN*2+1],
    	    left_eN_nDOF_trial_element = left_eN_global*nDOF_trial_element,
    	    right_eN_nDOF_trial_element = right_eN_global*nDOF_trial_element;
    	  double jac[nSpace*nSpace],
    	    jacDet,
    	    jacInv[nSpace*nSpace],
    	    boundaryJac[nSpace*(nSpace-1)],
    	    metricTensor[(nSpace-1)*(nSpace-1)],
    	    metricTensorDetSqrt,
    	    normal[3],
    	    x,y,z;
    	  //double G[nSpace*nSpace],G_dd_G,tr_G,h_phi,h_penalty;
	  
    	  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
    	    {
    	      ck.calculateMapping_elementBoundary(left_eN_global,
    						  left_ebN_element,
    						  kb,
    						  left_ebN_element*kb,
    						  mesh_dof,
    						  mesh_l2g,
    						  mesh_trial_trace_ref,
    						  mesh_grad_trial_trace_ref,
    						  boundaryJac_ref,
    						  jac,
    						  jacDet,
    						  jacInv,
    						  boundaryJac,
    						  metricTensor,
    						  metricTensorDetSqrt,
    						  normal_ref,
    						  normal,
    						  x,y,z);
    	      xArray_left[kb*3+0] = x;
    	      xArray_left[kb*3+1] = y;
    	      xArray_left[kb*3+2] = z;
    	      ck.calculateMapping_elementBoundary(right_eN_global,
    						  right_ebN_element,
    						  kb,
    						  right_ebN_element*kb,
    						  mesh_dof,
    						  mesh_l2g,
    						  mesh_trial_trace_ref,
    						  mesh_grad_trial_trace_ref,
    						  boundaryJac_ref,
    						  jac,
    						  jacDet,
    						  jacInv,
    						  boundaryJac,
    						  metricTensor,
    						  metricTensorDetSqrt,
    						  normal_ref,
    						  normal,
    						  x,y,z);
    	      xArray_right[kb*3+0] = x;
    	      xArray_right[kb*3+1] = y;
    	      xArray_right[kb*3+2] = z;
    	    }
    	  for  (int kb_left=0;kb_left<nQuadraturePoints_elementBoundary;kb_left++)
    	    {
    	      double errorNormMin = 1.0;
    	      for  (int kb_right=0;kb_right<nQuadraturePoints_elementBoundary;kb_right++)
    		{
    		  double errorNorm=0.0;
    		  for (int I=0;I<nSpace;I++)
    		    {
    		      errorNorm += fabs(xArray_left[kb_left*3+I]
    					-
    					xArray_right[kb_right*3+I]);
    		    }
    		  if (errorNorm < errorNormMin)
    		    {
    		      permutations[kb_right] = kb_left;
    		      errorNormMin = errorNorm;
    		    }
    		}
    	    }
    	  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
    	    {
    	      register int ebN_kb_nSpace = ebN*nQuadraturePoints_elementBoundary*nSpace+kb*nSpace;
    	      register double u_left=0.0,
    		v_left=0.0,
    		w_left=0.0,
    		u_right=0.0,
    		v_right=0.0,
    		w_right=0.0;
    	      register int left_kb = kb,
    		right_kb = permutations[kb],
    		left_ebN_element_kb_nDOF_test_element=left_ebN_element*left_kb*nDOF_test_element,
    		right_ebN_element_kb_nDOF_test_element=right_ebN_element*right_kb*nDOF_test_element;
    	      //
    	      //calculate the velocity solution at quadrature points on left and right
    	      //
    	      ck.valFromDOF(u_dof,&vel_l2g[left_eN_nDOF_trial_element],&vel_trial_trace_ref[left_ebN_element_kb_nDOF_test_element],u_left);
    	      ck.valFromDOF(v_dof,&vel_l2g[left_eN_nDOF_trial_element],&vel_trial_trace_ref[left_ebN_element_kb_nDOF_test_element],v_left);
    	      ck.valFromDOF(w_dof,&vel_l2g[left_eN_nDOF_trial_element],&vel_trial_trace_ref[left_ebN_element_kb_nDOF_test_element],w_left);
    	      //
    	      ck.valFromDOF(u_dof,&vel_l2g[right_eN_nDOF_trial_element],&vel_trial_trace_ref[right_ebN_element_kb_nDOF_test_element],u_right);
    	      ck.valFromDOF(v_dof,&vel_l2g[right_eN_nDOF_trial_element],&vel_trial_trace_ref[right_ebN_element_kb_nDOF_test_element],v_right);
    	      ck.valFromDOF(w_dof,&vel_l2g[right_eN_nDOF_trial_element],&vel_trial_trace_ref[right_ebN_element_kb_nDOF_test_element],w_right);
    	      //
    	      velocityAverage[ebN_kb_nSpace+0]=0.5*(u_left + u_right);
    	      velocityAverage[ebN_kb_nSpace+1]=0.5*(v_left + v_right);
    	      velocityAverage[ebN_kb_nSpace+2]=0.5*(w_left + w_right);
    	    }//ebNI
    	}
    }
  };//RBLES
  
  inline RBLES_base* newRBLES(int nSpaceIn,
		              int nQuadraturePoints_elementIn,
				int nDOF_mesh_trial_elementIn,
				int nDOF_trial_elementIn,
				int nDOF_test_elementIn,
				int nQuadraturePoints_elementBoundaryIn,
				int CompKernelFlag)
  {
    return proteus::chooseAndAllocateDiscretization<RBLES_base,RBLES,CompKernel>(nSpaceIn,
										   nQuadraturePoints_elementIn,
										   nDOF_mesh_trial_elementIn,
										   nDOF_trial_elementIn,
										   nDOF_test_elementIn,
										   nQuadraturePoints_elementBoundaryIn,
										   CompKernelFlag);
  }
}//proteus

#endif
