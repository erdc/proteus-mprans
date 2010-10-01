#include <iostream>

namespace Proteus
{
  template<class Model_Base, 
	   template<const int NSPACE, 
		    const int NDOF_MESH_TRIAL_ELEMENT, 
		    const int NDOF_TRIAL_ELEMENT, 
		    const int NDOF_TEST_ELEMENT> 
	   class ModelTemplate>
  Model_Base* chooseAndAllocateDiscretization(int nSpace, int nDOF_mesh, int nDOF_trial, int nDOF_test)
  {
    if (nSpace == 3)
      {
	if (nDOF_mesh == nDOF_trial)//iso-parametric
	  {
	    if (nDOF_mesh == 4)
	      return static_cast<Model_Base*>(new ModelTemplate<3,4,4,4>());
	    else if (nDOF_mesh == 8)
	      return static_cast<Model_Base*>(new ModelTemplate<3,8,8,8>());
	    else if (nDOF_mesh == 10)
	      return static_cast<Model_Base*>(new ModelTemplate<3,10,10,10>());
	    else if (nDOF_mesh == 27)
	      return static_cast<Model_Base*>(new ModelTemplate<3,27,27,27>());
	    else
	      abort();
	  }
	else if (nDOF_mesh == 4)//sub-parametric tets
	  {
	    if (nDOF_trial == 10)
	      return static_cast<Model_Base*>(new ModelTemplate<3,4,10,10>());
	    else
	      abort();
	  }
	else if (nDOF_mesh == 8)//sub-parametric hexes
	  {
	    if (nDOF_trial == 27)
	      return static_cast<Model_Base*>(new ModelTemplate<3,8,27,27>());
	    else
	      abort();
	  }
	else
	  abort();
      }
    else
      abort();
  }

  class MPRANS_Base
  {
  public:
    virtual void getRes(double* res, int nElements)=0;
  };

  template<const int NSPACE, const int NDOF_MESH_TRIAL_ELEMENT, const int NDOF_TRIAL_ELEMENT, const int NDOF_TEST_ELEMENT>
  class MPRANS : public MPRANS_Base
  {
  public:
    MPRANS()
    {
      std::cout<<"Hello World MPRANS()"<<std::endl;
    }
    virtual void getRes(double* res, int nElements)
    {
      std::cout<<"getRes "<<NSPACE<<" "<<NDOF_MESH_TRIAL_ELEMENT<<" "<<NDOF_TRIAL_ELEMENT<<" "<<NDOF_TEST_ELEMENT<<std::endl;
      for (int eN=0;eN<nElements;eN++)
	{
	  for (int i=0;i<NDOF_TRIAL_ELEMENT;i++)
	    for (int I=0;I<NSPACE;I++)
	      res[eN] = eN+i*3.0+I*4.0;
	}
    }
  };

  inline MPRANS_Base* newMPRANS(int nSpace, int nDOF_mesh, int nDOF_trial, int nDOF_test)
  {
    return chooseAndAllocateDiscretization<MPRANS_Base,MPRANS>(nSpace,nDOF_mesh,nDOF_trial,nDOF_test);
  }
}

