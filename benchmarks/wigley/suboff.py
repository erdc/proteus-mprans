#! /usr/bin/env python


def darpa2gen(npoints):
    """
    THIS PROGRAM CONTAINS FOLLOWING EQUATIONS FOR GENERATING
     OFFSETS IN FEET FOR DARPA2 MODEL WITH (FULL/MODEL) SCALE RATIO =24.

   INCLUDED ARE: 
     BOW EQ.	FOR	0.0 FT <- X <- 3.333333 FT, 
     PARALLEL MID-BODY EQ. FOR 3.333333 FT <- X <- 10.645833 FT, 
     AFTERBODY EQ.	FOR 10.645833 FT <- X <- 13.979167 FT, 
     AFTERBODY CAP EQ.	FOR 13-979167 FT <- X <- 14.291667 FT.
   AS SET UP HERE, OFFSETS ARE COMPUTED EVERY 0.1 FT. 
    (EXCEPT IN FIRST 0.5 FT, WHERE THEY ARE EVERY 6.01 FT)

    """
    from math import sqrt,pow
    import numpy
    x = numpy.zeros((npoints,),'d'); y = numpy.zeros((npoints,),'d')
    #constants
    rmax = 0.8333333
    xb   = 3.333333
    xm   = 10.645833
    xa   = 13.979167
    xc   = 14.291667
    cb1  = 1.126395191
    cb2  = 0.442874707
    cb3  = 1.0/2.1
    rh   = 0.1175
    k0   = 10.0
    k1   = 44.6244
    #
    xx = -0.01
    dx = 0.01
    i = 0
    while (i < npoints and xx < xc):
        np = i
        xx += dx
        if (xx >= 0.5): dx = 0.1
        if (xx >= xa): dx = 0.01
        if (xx < xb): #otherwise skip to 200
            #bow equation
            a = 0.3*xx - 1.0
            a3= a*a*a
            a4= a3*a
            b = 1.2*xx + 1.0
            r = cb1*xx*a4 + cb2*xx*xx*a3 + 1.0 - a4*b
            r = rmax*pow(r,cb3)
            x[i] = xx
            y[i] = r
        else: #goto 200
            if (xx < xm): #otherwise skip to 400
                #parallel mid-body equation
                x[i]=xx
                y[i]=rmax
            else: #goto 400
                if (xx < xa): #otherwise goto 600
                    #afterbody equation
                    xi = (13.979167 - xx)/3.333333
                    c1 = rh*rh
                    xipow = xi*xi
                    c2 = rh*k0*xipow                  
                    xipow *= xi #3
                    c3 = (20.0 - 20.0*rh*rh - 4.0*rh*k0 - 0.333333*k1)*xipow
                    xipow *= xi #4
                    c4 = (-45.0 + 45.0*rh*rh + 6.0*rh*k0 + k1)*xipow
                    xipow *= xi #5
                    c5 = (36.0  - 36.0*rh*rh - 4.0*rh*k0 - k1)*xipow
                    xipow *= xi #6
                    c6 = (-10.0 + 10.0*rh*rh +     rh*k0 + 0.333333*k1)*xipow
                    r  = rmax*sqrt((c1+c2+c3+c4+c5+c6))
                    x[i] = xx
                    y[i] = r
                else: #goto 600
                    if (xx < xc): #otherwise goto 1100
                        #afterbody cap equation
                        r = 1.0 - (3.2*xx - 44.733333)*(3.2*xx - 44.733333)
                        assert r >= 0.0, "negative square root in afterbody cap equation"
                        r = rh*rmax*sqrt(r)
                        x[i] = xx
                        y[i] = r
                    else:#1100
                        x[np] = xc
                        y[np] = 0.0
                    #end 1100 block
                #end 600 block
            #end 400 block
        #end 200 block
        i += 1
    #end loop 
    return x,y,np+1
#end darpagen2


if __name__ == '__main__':
    nx = 300
    x,y,np = darpa2gen(nx)

    fout = open('darpa2.dat','w')
    for i in range(np):
        fout.write("%12.5e %12.5e \n" % (x[i],y[i]))

    #
    fout.close()
    import matplotlib
    from matplotlib import pylab
    
    pylab.figure(1)
    pylab.plot(x[:np],y[:np])
    pylab.savefig('darpa2.png')
    pylab.show()
