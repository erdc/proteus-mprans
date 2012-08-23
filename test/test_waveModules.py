import waveModules_Matt as wm
import numpy as np

try:
    from matplotlib.pylab import plot,show           # for 2D plots (checking)            
    #from enthought.mayavi import mlab               # for 3D plots 
except:
    pass

""" Testing different analytic solution for wavetank benchmark case."""

#see nose docs for more complex testing

def test_Linear2D(showPlots=False):
    """ Testing the Linearized 2D interface (phi) propagation. """
    g = (0.0,0.0,-9.81)        # gravity
    L = (5.0,5.0,0.25)         # tank dimensions

    # Water                                                                  
    rho_0 = 998.2
    nu_0  = 1.004e-6

    # Air                                                                    
    rho_1 = 1.205
    nu_1  = 1.500e-5
 
    inflowHeightMean = 0.5*L[2]
    inflowVelocityMean = (0.2,0.0,0.0)
    waveLength = 5*inflowHeightMean

    A = 0.1                # amplitude
    k = 2*np.pi/waveLength
    h = L[2]
    omega = np.sqrt(-g[2]*k*np.tanh(k*h))
    period = 2*np.pi/omega

    # Collocation point
    N = 100
    x = np.linspace(0,L[0],N)
    t = np.linspace(0,period,N)
    [xx, tt] = np.meshgrid(x,t) 
    
    # Wave Field Object
    # --> waveTest = wm.Linear2D(A,omega,k,h,rho_0,rho_1)
    waveTest = wm.WaveGroup(A,omega,k,h,rho_0,rho_1)
    # --> waveTest = wm.Solitary(A,omega,k,h,rho_0,rho_1)

    result = waveTest.height(x,t[90])
    #correctResult = np.real( A*np.exp(1j*(k*x - omega*t[1])) )

    # Plot result if appropriate
    #assert correctResult.all() == result.all(), "Linear2D.height returned %f should be %f" % (result,correctResult)

    try:
        plot(x, result, 'b', linewidth=2)
        if showPlots:
            show()
    except:
        pass


if __name__ == '__main__':
    print "The program name is: ", __name__
    test_Linear2D(showPlots=True)
    #-->test_WaveGroup()
    #-->test_Solitary()
