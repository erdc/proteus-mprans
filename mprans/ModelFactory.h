#ifndef MODELFACTORY_H
#define MODELFACTORY_H

namespace proteus
{
  template<class Model_Base, 
	   template<class CompKernelType,
		    int nSpace,
		    int nQuadraturePoints_element,
		    int nDOF_mesh_trial_element,
		    int nDOF_trial_element,
		    int nDOF_test_element,
		    int nQuadraturePoints_elementBoundary>
	   class ModelTemplate,
	   template<int nSpace,
		    int nDOF_mesh_trial_element,
		    int nDOF_trial_element,
		    int nDOF_test_element>
	   class CompKernelTemplate>
  Model_Base* chooseAndAllocateDiscretization(int nSpaceIn,
					      int nQuadraturePoints_elementIn,
					      int nDOF_mesh_trial_elementIn,
					      int nDOF_trial_elementIn,
					      int nDOF_test_elementIn,
					      int nQuadraturePoints_elementBoundaryIn,
					      int CompKernelFlag)//0=Parametric
	   {
	     if (CompKernelFlag == 0)
	       {
		 if (nSpaceIn == 3)
		   {
		     if (nDOF_mesh_trial_elementIn == nDOF_trial_elementIn)//iso-parametric
		       {
			 if (nDOF_mesh_trial_elementIn == 4)
			   {
			     if (nQuadraturePoints_elementIn == 5)
			       return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,4,4,4>,3,5,4,4,4,4>());
			     else
			       abort();
			   }
			 else if (nDOF_mesh_trial_elementIn == 8)
			   {
			     abort();//ido todo, add q1 iso
			   }
			 else if (nDOF_mesh_trial_elementIn == 10)
			   {
			     if (nQuadraturePoints_elementIn == 15)
			       return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,10,10,10>,3,15,10,10,10,7>());
			     else
			       abort();
			   }
			 else if (nDOF_mesh_trial_elementIn == 27)
			   {
			     abort();//ido todo, add q2 iso
			   }
			 else
			   abort();
		       }
		     else if (nDOF_mesh_trial_elementIn == 4)//sub-parametric tets
		       {
			 if (nDOF_trial_elementIn == 10)
			   return static_cast<Model_Base*>(new ModelTemplate<CompKernelTemplate<3,4,10,10>,3,15,4,10,10,7>());
			 else
			   abort();
		       }
		     else if (nDOF_mesh_trial_elementIn == 8)//sub-parametric hexes
		       {
			 abort();//todo ito add q2 on hex mesh
		       }
		     else
		       abort();
		   }
		 else
		   abort();
	       }
	     else
	       abort();
	   }
}
#endif
