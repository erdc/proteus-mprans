#include "RBLES.h"
#include <iostream>
#include <cassert>

extern "C" void calculateResidual_RBLES(//element
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
					   int nElements_global,
                                           double useRBLES,
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
				  int* p_IBC, int* u_IBC, int* v_IBC, int* w_IBC,
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
					   double* flux)
{
  CompKernel<nSpace,nDOF_mesh_trial_element,nDOF_trial_element,nDOF_test_element> ck;

//  for(int i=0;i<120;i++)
//	  std::cout<<u_IBC[i]<<std::endl;
  //
  //loop over elements to compute volume integrals and load them into element and global residual
  //
  double globalConservationError=0.0;
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
	    rho,mu,H_rho,d_rho,H_mu,d_mu,norm_n;
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
 h_phi = 0.1;
	  eps_rho = epsFact_rho*h_phi;
	  eps_mu = epsFact_mu*h_phi;
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

	
          //calculate tau and tau*Res
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




          double Dt=(dmom_u_acc_u_t/dmom_u_acc_u);
          double v_d_Gv= 0.0;
	  double vel[nSpace] = {u,v,w};
          for(int I=0;I<nSpace;I++)
            for (int J=0;J<nSpace;J++)
	      v_d_Gv += vel[I]*G[I*nSpace+J]*vel[J];
  
          tau_0 = 1.0/sqrt(Ct_sge*rho*rho*Dt*Dt + rho*rho*v_d_Gv + Cd_sge*mu*mu*G_dd_G);
          tau_1 = 1.0/(tr_G*tau_0);
          //tau_0 = 0.0; 
	  //tau_1 = 0.0;
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
	  mu += C_dc*norm_Rv/sqrt(G_dd_G);

          // 
          //update element residual 
          // 
          for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      register int i_nSpace=i*nSpace;	



/*
              elementResidual_u[i] += (u - g[0])*vel_test_dV[i];	
              elementResidual_v[i] += (v - g[1])*vel_test_dV[i];	
              elementResidual_w[i] += (w - g[2])*vel_test_dV[i];					    


              elementResidual_p[i] += grad_p[0]*p_grad_test_dV[i_nSpace+0]
	                            + grad_p[1]*p_grad_test_dV[i_nSpace+1]
	                            + grad_p[2]*p_grad_test_dV[i_nSpace+2] + (p+1.0)*p_test_dV[i];*/
				    

	      elementResidual_p[i] += -pdeResidual_p*p_test_dV[i]
	                            + subgrid_u*p_grad_test_dV[i_nSpace+0]
	                            + subgrid_v*p_grad_test_dV[i_nSpace+1]
	                            + subgrid_w*p_grad_test_dV[i_nSpace+2];

	      elementResidual_u[i] += rho *(mom_u_acc_t+ u*grad_u[0] + v*grad_u[1] + w*grad_u[2] - g[0])*vel_test_dV[i]				        
			               + mu*(grad_u[0] + grad_u[0])*vel_grad_test_dV[i_nSpace+0]
			               + mu*(grad_u[1] + grad_v[0])*vel_grad_test_dV[i_nSpace+1]				   
			               + mu*(grad_u[2] + grad_w[0])*vel_grad_test_dV[i_nSpace+2]	
	                               -(p+subgrid_p)*vel_grad_test_dV[i_nSpace+0]     	 
		                        -rho*vel_grad_test_dV[i_nSpace+0] * (subgrid_u*u) //+ u*subgrid_u + subgrid_u*subgrid_u)
		                        -rho*vel_grad_test_dV[i_nSpace+1] * (subgrid_u*v) //+ u*subgrid_v + subgrid_u*subgrid_v)
		                        -rho*vel_grad_test_dV[i_nSpace+2] * (subgrid_u*w);// + u*subgrid_w + subgrid_u*subgrid_w);
		            
	      elementResidual_v[i] += rho *(mom_v_acc_t+ u*grad_v[0] + v*grad_v[1] + w*grad_v[2] - g[1])*vel_test_dV[i]
			               + mu*(grad_v[0] + grad_u[1])*vel_grad_test_dV[i_nSpace+0]
			               + mu*(grad_v[1] + grad_v[1])*vel_grad_test_dV[i_nSpace+1]				   
			               + mu*(grad_v[2] + grad_w[1])*vel_grad_test_dV[i_nSpace+2]	
	                               -(p+subgrid_p)*vel_grad_test_dV[i_nSpace+1] 	                                  
		                       -rho*vel_grad_test_dV[i_nSpace+0] * (subgrid_v*u) //+ v*subgrid_u + subgrid_v*subgrid_u)
		                       -rho*vel_grad_test_dV[i_nSpace+1] * (subgrid_v*v) //+ v*subgrid_v + subgrid_v*subgrid_v)
		                       -rho*vel_grad_test_dV[i_nSpace+2] * (subgrid_v*w);// + v*subgrid_w + subgrid_v*subgrid_w);


	      elementResidual_w[i] += rho *(mom_w_acc_t + u*grad_w[0] + v*grad_w[1] + w*grad_w[2] - g[2])*vel_test_dV[i]
			               + mu*(grad_w[0] + grad_u[2])*vel_grad_test_dV[i_nSpace+0]
			               + mu*(grad_w[1] + grad_v[2])*vel_grad_test_dV[i_nSpace+1]				   
			               + mu*(grad_w[2] + grad_w[2])*vel_grad_test_dV[i_nSpace+2]		 
	                               -(p+subgrid_p)*vel_grad_test_dV[i_nSpace+2] 
		                        -rho*vel_grad_test_dV[i_nSpace+0] * (subgrid_w*u) //+ w*subgrid_u + subgrid_w*subgrid_u)
		                        -rho*vel_grad_test_dV[i_nSpace+1] * (subgrid_w*v) //+ w*subgrid_v + subgrid_w*subgrid_v)
		                        -rho*vel_grad_test_dV[i_nSpace+2] * (subgrid_w*w);// + w*subgrid_w + subgrid_w*subgrid_w);
           }//i
	}
      // 
      //load element into global residual and save element residual
      //
      for(int i=0;i<nDOF_test_element;i++) 
        { 
          register int eN_i=eN*nDOF_test_element+i;

          if (p_IBC[p_l2g[eN_i]]   == 1)  elementResidual_p[i]=0.0;
          if (u_IBC[vel_l2g[eN_i]] == 1)  elementResidual_u[i]=0.0;
          if (v_IBC[vel_l2g[eN_i]] == 1)  elementResidual_v[i]=0.0;
          if (w_IBC[vel_l2g[eN_i]] == 1)  elementResidual_w[i]=0.0;


	  
          q_elementResidual_p[eN_i]+=elementResidual_p[i];
          q_elementResidual_u[eN_i]+=elementResidual_u[i];
          q_elementResidual_v[eN_i]+=elementResidual_v[i];
          q_elementResidual_w[eN_i]+=elementResidual_w[i];
          globalResidual[offset_p+stride_p*p_l2g[eN_i]]  += elementResidual_p[i];
          globalResidual[offset_u+stride_u*vel_l2g[eN_i]]+= elementResidual_u[i];
          globalResidual[offset_v+stride_v*vel_l2g[eN_i]]+= elementResidual_v[i];
          globalResidual[offset_w+stride_w*vel_l2g[eN_i]]+= elementResidual_w[i];
        }//i

    }//elements

}

extern "C" void calculateJacobian_RBLES(//element
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
					   int nElements_global,
                                           double useRBLES,
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
				  int* p_IBC, int* u_IBC, int* v_IBC, int* w_IBC,
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
  CompKernel<nSpace,nDOF_mesh_trial_element,nDOF_trial_element,nDOF_test_element> ck;

//  for(int i=0;i<120;i++)
//	  std::cout<<u_IBC[i]<<std::endl;


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
		 
	 h_phi =0.1;
	  eps_rho = epsFact_rho*h_phi;
	  eps_mu = epsFact_mu*h_phi;
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
	                             + u*vel_grad_trial[j_nSpace+0] 
	                             + v*vel_grad_trial[j_nSpace+1] 
	                             + w*vel_grad_trial[j_nSpace+2]);

	      dpdeResidual_v_p[j]=p_grad_trial[j_nSpace+1];
	      dpdeResidual_v_v[j]=rho*(dmom_v_acc_v_t*vel_trial_ref[j]
	                             + u*vel_grad_trial[j_nSpace+0] 
	                             + v*vel_grad_trial[j_nSpace+1] 
	                             + w*vel_grad_trial[j_nSpace+2]);
	                             
	      dpdeResidual_w_p[j]=p_grad_trial[j_nSpace+2];
	      dpdeResidual_w_w[j]=rho*(dmom_w_acc_w_t*vel_trial_ref[j]
	                             + u*vel_grad_trial[j_nSpace+0] 
	                             + v*vel_grad_trial[j_nSpace+1] 
	                             + w*vel_grad_trial[j_nSpace+2]);
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
	  mu += C_dc*norm_Rv/sqrt(G_dd_G);


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
		  elementJacobian_u_u[i][j] += dpdeResidual_u_u[j]*vel_test_dV[i];	 
		                        -rho*vel_grad_test_dV[i_nSpace+0] * (dsubgrid_u_u[j]*u)// + vel_trial_ref[j]*subgrid_u)
		                        -rho*vel_grad_test_dV[i_nSpace+1] * (dsubgrid_u_u[j]*v)// + vel_trial_ref[j]*subgrid_v)
		                        -rho*vel_grad_test_dV[i_nSpace+2] * (dsubgrid_u_u[j]*w);// + vel_trial_ref[j]*subgrid_w);

		  elementJacobian_u_p[i][j] += -p_trial_ref[j]*vel_grad_test_dV[i_nSpace+0];	                                  
		                        -rho*vel_grad_test_dV[i_nSpace+0] * (dsubgrid_u_p[j]*u)// + v*dsubgrid_u_p[j])
		                        -rho*vel_grad_test_dV[i_nSpace+1] * (dsubgrid_u_p[j]*v)//+ v*dsubgrid_v_p[j])
		                        -rho*vel_grad_test_dV[i_nSpace+2] * (dsubgrid_u_p[j]*w);// + v*dsubgrid_w_p[j]);

		  elementJacobian_u_u[i][j] += -dsubgrid_p_u[j]*vel_grad_test_dV[i_nSpace+0]; 		  
		  elementJacobian_u_v[i][j] += -dsubgrid_p_v[j]*vel_grad_test_dV[i_nSpace+0];
		  elementJacobian_u_w[i][j] += -dsubgrid_p_w[j]*vel_grad_test_dV[i_nSpace+0];		

		  elementJacobian_u_u[i][j] += 2.0*mu*vel_grad_trial[j_nSpace+0]*vel_grad_test_dV[i_nSpace+0]
		                                +  mu*vel_grad_trial[j_nSpace+1]*vel_grad_test_dV[i_nSpace+1]			   
		                                +  mu*vel_grad_trial[j_nSpace+2]*vel_grad_test_dV[i_nSpace+2];

		  elementJacobian_u_v[i][j] += mu*vel_grad_trial[j_nSpace+0]*vel_grad_test_dV[i_nSpace+1];				   
		  elementJacobian_u_w[i][j] += mu*vel_grad_trial[j_nSpace+0]*vel_grad_test_dV[i_nSpace+2];


		  


          // v 
		  elementJacobian_v_v[i][j] += dpdeResidual_v_v[j]*vel_test_dV[i];		 
		                        -rho*vel_grad_test_dV[i_nSpace+0] * (dsubgrid_v_v[j]*u)// + vel_trial_ref[j]*subgrid_u)
		                        -rho*vel_grad_test_dV[i_nSpace+1] * (dsubgrid_v_v[j]*v)// + vel_trial_ref[j]*subgrid_v)
		                        -rho*vel_grad_test_dV[i_nSpace+2] * (dsubgrid_v_v[j]*w);// + vel_trial_ref[j]*subgrid_w);
					
		  elementJacobian_v_p[i][j] += -p_trial_ref[j]*vel_grad_test_dV[i_nSpace+1]	;                                  
		                        -rho*vel_grad_test_dV[i_nSpace+0] * (dsubgrid_v_p[j]*u)// + v*dsubgrid_u_p[j])
		                        -rho*vel_grad_test_dV[i_nSpace+1] * (dsubgrid_v_p[j]*v)// + v*dsubgrid_v_p[j])
		                        -rho*vel_grad_test_dV[i_nSpace+2] * (dsubgrid_v_p[j]*w);// + v*dsubgrid_w_p[j]);

		  elementJacobian_v_u[i][j] += -dsubgrid_p_u[j]*vel_grad_test_dV[i_nSpace+1]; 		  
		  elementJacobian_v_v[i][j] += -dsubgrid_p_v[j]*vel_grad_test_dV[i_nSpace+1];
		  elementJacobian_v_w[i][j] += -dsubgrid_p_w[j]*vel_grad_test_dV[i_nSpace+1];		  

		  elementJacobian_v_u[i][j] += mu*vel_grad_trial[j_nSpace+1]*vel_grad_test_dV[i_nSpace+0];
		  elementJacobian_v_v[i][j] +=     mu*vel_grad_trial[j_nSpace+0]*vel_grad_test_dV[i_nSpace+0]
		                             + 2.0*mu*vel_grad_trial[j_nSpace+1]*vel_grad_test_dV[i_nSpace+1]				   
		                                +  mu*vel_grad_trial[j_nSpace+2]*vel_grad_test_dV[i_nSpace+2];				   
		  elementJacobian_v_w[i][j] += mu*vel_grad_trial[j_nSpace+1]*vel_grad_test_dV[i_nSpace+2];
		  		  


		 
          // w 
		  elementJacobian_w_w[i][j] += dpdeResidual_w_w[j]*vel_test_dV[i];		 
		                        -rho*vel_grad_test_dV[i_nSpace+0] * (dsubgrid_w_w[j]*u)// + vel_trial_ref[j]*subgrid_u)
		                        -rho*vel_grad_test_dV[i_nSpace+1] * (dsubgrid_w_w[j]*v)// + vel_trial_ref[j]*subgrid_v)
		                        -rho*vel_grad_test_dV[i_nSpace+2] * (dsubgrid_w_w[j]*w);// + vel_trial_ref[j]*subgrid_w);
		  
		  elementJacobian_w_p[i][j] += -p_trial_ref[j]*vel_grad_test_dV[i_nSpace+2];	                                  
		                        -rho*vel_grad_test_dV[i_nSpace+0] * (dsubgrid_w_p[j]*u)// + v*dsubgrid_u_p[j])
		                        -rho*vel_grad_test_dV[i_nSpace+1] * (dsubgrid_w_p[j]*v)// + v*dsubgrid_v_p[j])
		                        -rho*vel_grad_test_dV[i_nSpace+2] * (dsubgrid_w_p[j]*w);// + v*dsubgrid_w_p[j]);

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
	  
	  
	  
          if (p_IBC[p_l2g[eN_i]]   == 1) {
	    for (int j=0;j<nDOF_trial_element;j++)
	    {
	      elementJacobian_p_p[i][j] = 0.0;
	      elementJacobian_p_u[i][j] = 0.0;
	      elementJacobian_p_v[i][j] = 0.0;
	      elementJacobian_p_w[i][j] = 0.0;

	      elementJacobian_p_p[j][i] = 0.0;
	      elementJacobian_u_p[j][i] = 0.0;
	      elementJacobian_v_p[j][i] = 0.0;
	      elementJacobian_w_p[j][i] = 0.0;	      	      	      

	    }
	    elementJacobian_p_p[i][i] = 1.0;
	  }
	  
          if (u_IBC[vel_l2g[eN_i]] == 1)	  {
	    for (int j=0;j<nDOF_trial_element;j++)
	    {
	      elementJacobian_u_p[i][j] = 0.0;
	      elementJacobian_u_u[i][j] = 0.0;
	      elementJacobian_u_v[i][j] = 0.0;
	      elementJacobian_u_w[i][j] = 0.0;

	      elementJacobian_p_u[j][i] = 0.0;
	      elementJacobian_u_u[j][i] = 0.0;
	      elementJacobian_v_u[j][i] = 0.0;
	      elementJacobian_w_u[j][i] = 0.0;	      	      	      

	    }
	    elementJacobian_u_u[i][i] = 1.0;
	  }
          if (v_IBC[vel_l2g[eN_i]] == 1)	  {
	    for (int j=0;j<nDOF_trial_element;j++)
	    {
	      elementJacobian_v_p[i][j] = 0.0;
	      elementJacobian_v_u[i][j] = 0.0;
	      elementJacobian_v_v[i][j] = 0.0;
	      elementJacobian_v_w[i][j] = 0.0;

	      elementJacobian_p_v[j][i] = 0.0;
	      elementJacobian_u_v[j][i] = 0.0;
	      elementJacobian_v_v[j][i] = 0.0;
	      elementJacobian_w_v[j][i] = 0.0;	      	      	      

	    }
	    elementJacobian_v_v[i][i] = 1.0;
	  }
          if (w_IBC[vel_l2g[eN_i]] == 1)	  {
	    for (int j=0;j<nDOF_trial_element;j++)
	    {
	      elementJacobian_w_p[i][j] = 0.0;
	      elementJacobian_w_u[i][j] = 0.0;
	      elementJacobian_w_v[i][j] = 0.0;
	      elementJacobian_w_w[i][j] = 0.0;

	      elementJacobian_p_w[j][i] = 0.0;
	      elementJacobian_u_w[j][i] = 0.0;
	      elementJacobian_v_w[j][i] = 0.0;
	      elementJacobian_w_w[j][i] = 0.0;	      	      	      

	    }
	    elementJacobian_w_w[i][i] = 1.0;
	  }

        }


   /*   for (int i=0;i<nDOF_test_element;i++)
	{

	  register int eN_i = eN*nDOF_test_element+i;
	  
          if (p_IBC[p_l2g[eN_i]]   == 1) {
	    for (int j=0;j<nDOF_trial_element;j++)
	    {
	      elementJacobian_p_p[j][i] = 0.0;
	      elementJacobian_p_u[j][i] = 0.0;
	      elementJacobian_p_v[j][i] = 0.0;
	      elementJacobian_p_w[j][i] = 0.0;

	      elementJacobian_p_p[i][j] = 0.0;
	      elementJacobian_u_p[i][j] = 0.0;
	      elementJacobian_v_p[i][j] = 0.0;
	      elementJacobian_w_p[i][j] = 0.0;	      	      	      

	    }
	    elementJacobian_p_p[i][i] = 1.0;
	  }
          if (u_IBC[vel_l2g[eN_i]] == 1)	  {
	    for (int j=0;j<nDOF_trial_element;j++)
	    {
	      elementJacobian_u_p[j][i] = 0.0;
	      elementJacobian_u_u[j][i] = 0.0;
	      elementJacobian_u_v[j][i] = 0.0;
	      elementJacobian_u_w[j][i] = 0.0;

	      elementJacobian_p_u[i][j] = 0.0;
	      elementJacobian_u_u[i][j] = 0.0;
	      elementJacobian_v_u[i][j] = 0.0;
	      elementJacobian_w_u[i][j] = 0.0;	      	      	      

	    }
	    elementJacobian_u_u[i][i] = 1.0;
	  }
          if (v_IBC[vel_l2g[eN_i]] == 1)	  {
	    for (int j=0;j<nDOF_trial_element;j++)
	    {
	      elementJacobian_v_p[j][i] = 0.0;
	      elementJacobian_v_u[j][i] = 0.0;
	      elementJacobian_v_v[j][i] = 0.0;
	      elementJacobian_v_w[j][i] = 0.0;

	      elementJacobian_p_v[i][j] = 0.0;
	      elementJacobian_u_v[i][j] = 0.0;
	      elementJacobian_v_v[i][j] = 0.0;
	      elementJacobian_w_v[i][j] = 0.0;	      	      	      

	    }
	    elementJacobian_v_v[i][i] = 1.0;
	  }
          if (w_IBC[vel_l2g[eN_i]] == 1)	  {
	    for (int j=0;j<nDOF_trial_element;j++)
	    {
	      elementJacobian_w_p[j][i] = 0.0;
	      elementJacobian_w_u[j][i] = 0.0;
	      elementJacobian_w_v[j][i] = 0.0;
	      elementJacobian_w_w[j][i] = 0.0;

	      elementJacobian_p_w[i][j] = 0.0;
	      elementJacobian_u_w[i][j] = 0.0;
	      elementJacobian_v_w[i][j] = 0.0;
	      elementJacobian_w_w[i][j] = 0.0;	      	      	      

	    }
	    elementJacobian_w_w[i][i] = 1.0;
	  }

        }*/


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
      // //
      // //debug element jacobian
      // //
      // std::cout<<"element jacobian"<<std::endl;
       /*for (int i=0;i<nDOF_test_element;i++)
       	{
       	  for (int j=0;j<nDOF_trial_element;j++)
       	    {
       	      std::cout <<" ---------------------------------------- "<<std::endl;
	      
	      std::cout << elementJacobian_p_p[i][j]<<"  ";
       	      std::cout << elementJacobian_p_u[i][j]<<"  ";
       	      std::cout << elementJacobian_p_v[i][j]<<"  ";
       	      std::cout << elementJacobian_p_w[i][j]<<std::endl;

       	      std::cout << elementJacobian_u_p[i][j]<<"  ";
      	      std::cout << elementJacobian_u_u[i][j]<<"  ";
       	      std::cout << elementJacobian_u_v[i][j]<<"  ";
       	      std::cout << elementJacobian_u_w[i][j]<<std::endl;

       	      std::cout << elementJacobian_v_p[i][j]<<"  ";
       	      std::cout << elementJacobian_v_u[i][j]<<"  ";
       	      std::cout << elementJacobian_v_v[i][j]<<"  ";
       	      std::cout << elementJacobian_v_w[i][j]<<std::endl;

       	      std::cout << elementJacobian_w_p[i][j]<<"  ";
       	      std::cout << elementJacobian_w_u[i][j]<<"  ";
       	      std::cout << elementJacobian_w_v[i][j]<<"  ";
       	      std::cout << elementJacobian_w_w[i][j]<<std::endl;
       	    }//j
       	}//i*/
    }//elements

}//computeJacobian

extern "C" void calculateVelocityAverage_RBLES(int* permutations,
						  int nExteriorElementBoundaries_global,
						  int* exteriorElementBoundariesArray,
						  int nInteriorElementBoundaries_global,
						  int* interiorElementBoundariesArray,
						  int* elementBoundaryElementsArray,
						  int* elementBoundaryLocalElementBoundariesArray,
						  int* vel_l2g, 
						  double* u_dof, double* v_dof, double* w_dof,
						  double* vel_trial_ref,
						  double* ebqe_velocity,
						  double* velocityAverage)
{
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
	right_ebN_element = elementBoundaryLocalElementBoundariesArray[ebN*2+1];

      for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	{ 
	  register int ebN_kb_nSpace = ebN*nQuadraturePoints_elementBoundary*nSpace+kb*nSpace;
	  register double u_left=0.0,
	    v_left=0.0,
	    w_left=0.0,
	    u_right=0.0,
	    v_right=0.0,
	    w_right=0.0;
	  register int left_kb = permutations[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary+
					      left_ebN_element*nQuadraturePoints_elementBoundary+
					      kb],
	    right_kb = permutations[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary+
				    right_ebN_element*nQuadraturePoints_elementBoundary+
				    kb];
	  // 
	  //calculate the velocity solution at quadrature points on left and right
	  // 
	  for (int j=0;j<nDOF_trial_element;j++) 
	    { 
	      int left_eN_j = left_eN_global*nDOF_trial_element+j;
	      int left_ebN_kb_j = left_ebN_element*nQuadraturePoints_elementBoundary*nDOF_trial_element + 
		left_kb*nDOF_trial_element +
		j;
	      int right_eN_j = right_eN_global*nDOF_trial_element+j;
	      int right_ebN_kb_j = right_ebN_element*nQuadraturePoints_elementBoundary*nDOF_trial_element + 
		right_kb*nDOF_trial_element +
		j;
	      u_left += u_dof[vel_l2g[left_eN_j]]*vel_trial_ref[left_ebN_kb_j]; 
	      v_left += v_dof[vel_l2g[left_eN_j]]*vel_trial_ref[left_ebN_kb_j]; 
	      w_left += w_dof[vel_l2g[left_eN_j]]*vel_trial_ref[left_ebN_kb_j]; 
	      u_right += u_dof[vel_l2g[right_eN_j]]*vel_trial_ref[right_ebN_kb_j]; 
	      v_right += v_dof[vel_l2g[right_eN_j]]*vel_trial_ref[right_ebN_kb_j]; 
	      w_right += w_dof[vel_l2g[right_eN_j]]*vel_trial_ref[right_ebN_kb_j]; 
	    }
	  velocityAverage[ebN_kb_nSpace+0]=0.5*(u_left + u_right);
	  velocityAverage[ebN_kb_nSpace+1]=0.5*(v_left + v_right);
	  velocityAverage[ebN_kb_nSpace+2]=0.5*(w_left + w_right);
	}//ebNI
    }
}
