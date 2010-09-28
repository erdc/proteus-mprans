import mprTest
import numpy

if __name__ == '__main__':
    res = numpy.zeros((100000,),"d")
    p1 = mprTest.PyMPR(3,4,4,4)#tets
    p1.getRes(res)
    print "iso linear tet",res
    p1 = mprTest.PyMPR(3,10,10,10)#tets
    p1.getRes(res)
    print "iso quadratic tet",res
    p1 = mprTest.PyMPR(3,4,10,10)#tets
    p1.getRes(res)
    print "quadratic on linear tet",res
    q1 = mprTest.PyMPR(3,8,8,8)#hex
    q1.getRes(res)
    print "iso bilinear hex",res
    q1 = mprTest.PyMPR(3,27,27,27)#hex
    q1.getRes(res)
    print "iso biquadratic hex",res
    q1 = mprTest.PyMPR(3,8,27,27)#hex
    q1.getRes(res)
    print "biquadratic on bilinear hex",res
