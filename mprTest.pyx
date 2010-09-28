import numpy
cimport numpy

cdef extern from "TestTemplate.h" namespace "Proteus":
   cdef cppclass MPRANS_Base:
      void getRes(double* res,int nElements)
   cdef cppclass MPRANS[N1,N2,N3,N4]:
      MPRANS()
      void getRes(double* res,int nElements)

cdef extern from *:
  ctypedef int N1 "1"
  ctypedef int N2 "2"
  ctypedef int N3 "3"
  ctypedef int N4 "4"
  ctypedef int N5 "5"
  ctypedef int N6 "6"
  ctypedef int N7 "7"
  ctypedef int N8 "8"
  ctypedef int N9 "9"
  ctypedef int N10 "10"
  ctypedef int N27 "27"

cdef class PyMPR:
   cdef MPRANS_Base* thisptr
   def __cinit__(self,nSpace,nDOF_mesh,nDOF_trial,nDOF_test):
       if nSpace == 3:
           if nDOF_mesh == nDOF_trial:#iso-parametric
               if nDOF_mesh == 4:
                   self.thisptr = <MPRANS_Base*> new MPRANS[N3,N4,N4,N4]()
               elif nDOF_mesh == 8:
                   self.thisptr = <MPRANS_Base*> new MPRANS[N3,N8,N8,N8]()
               elif nDOF_mesh == 10:
                   self.thisptr = <MPRANS_Base*> new MPRANS[N3,N10,N10,N10]()
               elif nDOF_mesh == 27:
                   self.thisptr = <MPRANS_Base*> new MPRANS[N3,N27,N27,N27]()
               else:
                   raise RuntimeError
           elif nDOF_mesh == 4:#sub-parametric tets
               if nDOF_trial == 10:
                   self.thisptr = <MPRANS_Base*> new MPRANS[N3,N4,N10,N10]()
               else:
                   raise RuntimeError
           elif nDOF_mesh == 8:#sub-parametric hexes
               if nDOF_trial == 27:
                   self.thisptr = <MPRANS_Base*> new MPRANS[N3,N8,N27,N27]()
               else:
                   raise RuntimeError
           else:
               raise RuntimeError
       else:
           raise RuntimeError

   def __dealloc__(self):
       del self.thisptr
   def getRes(self,numpy.ndarray[numpy.double_t,ndim=1] res):
       self.thisptr.getRes(<double*>res.data,res.shape[0])
