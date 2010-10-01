import numpy
cimport numpy

cdef extern from "TestTemplate.h" namespace "Proteus":
   cdef cppclass MPRANS_Base:
      void getRes(double* res,
                  int nElements)
   MPRANS_Base* newMPRANS(int nSpace, 
                          int nDOF_mesh, 
                          int nDOF_trial, 
                          int nDOF_test)

cdef class PyMPR:
   cdef MPRANS_Base* thisptr
   def __cinit__(self,nSpace,nDOF_mesh,nDOF_trial,nDOF_test):
       self.thisptr = newMPRANS(nSpace,nDOF_mesh,nDOF_trial,nDOF_test)
   def __dealloc__(self):
       del self.thisptr
   def getRes(self,numpy.ndarray[numpy.double_t,ndim=1] res):
       self.thisptr.getRes(<double*>res.data,res.shape[0])
