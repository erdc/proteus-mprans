#include <iostream>

namespace Proteus
{
  class MPRANS_Base
  {
  public:
    virtual void getRes(double* res, int nElements)=0;
  };
  template<const int NSPACE, const int NDOF_MESH_TRIAL_ELEMENT, const int NDOF_TRIAL_ELEMENT, const int NDOF_TEST_ELEMENT>
  class MPRANS
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
}
