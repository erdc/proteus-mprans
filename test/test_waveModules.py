import waveModules_Matt as wm
import JONSWAP_p as JS
import numpy as np

try:
    from matplotlib.pylab import plot,show           # for 2D plots (checking)            
    #from enthought.mayavi import mlab               # for 3D plots 
    import animation # in ~/proteus-mprans/.../wavetank/ as my version of matplotlib doesn't have it
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import mpl_toolkits.mplot3d.axes3d as p3
except:
    pass

""" Testing different analytic solution for wavetank benchmark case."""


g = (0.0,0.0,-9.81)        # gravity
L = (20.0,0.25,1.0)         # tank dimensions

# Water                                                                  
rho_0 = 998.2
nu_0  = 1.004e-6

# Air                                                                    
rho_1 = 1.205
nu_1  = 1.500e-5
 
inflowHeightMean = 0.5*L[2]
inflowVelocityMean = (0.0,0.0,0.0)
waveLength = 5*inflowHeightMean
amplitude = 0.1*inflowHeightMean

#see nose docs for more complex testing

def test_Linear2D(showPlots=False):
    """ Testing the Linearized 2D interface (phi) propagation. """
    A = amplitude                # amplitude
    k = (2*np.pi/waveLength,0.0,0.0)
    h = L[2]
    omega = np.sqrt(-g[2]*k[0]*np.tanh(k[0]*h))
    period = 2*np.pi/omega

    # Collocation point
    N = 100
    x = [np.linspace(0,L[0],N), 0.0, 0.0]
    t = np.linspace(0,period,N)
    [xx, tt] = np.meshgrid(x,t) 
    
    # Wave Field Object
    waveTest = wm.Linear2D(A,omega,k,h,rho_0,rho_1)

    result = waveTest.height(x,t[90])
    #correctResult = np.real( A*np.exp(1j*(np.inner(k,x)- omega*t[1])) )

    # Plot result if appropriate
    #assert correctResult.all() == result.all(), "Linear2D.height returned %f should be %f" % (result,correctResult)

    try:
        plot(x[0], result, 'b', linewidth=2)
        if showPlots:
            show()
    except:
        pass


def test_WaveGroup(showPlots=False):
    """ Testing the Linearized 2D interface (phi) propagation. """
    A = amplitude                # amplitude
    k = (2*np.pi/waveLength,0.0,0.0)
    h = L[2]
    omega = np.sqrt(-g[2]*k[0]*np.tanh(k[0]*h))
    period = 2*np.pi/omega

    # Collocation point
    N = 100
    x = [np.linspace(0,L[0],N), 0.0, 0.0]
    t = np.linspace(0,period,N)
    [xx, tt] = np.meshgrid(x,t) 
    
    # Wave Field Object
    waveTest = wm.WaveGroup(A,omega,k,h,rho_0,rho_1)

    result = waveTest.height(x,t[90])
    #correctResult = np.real( A*np.exp(1j*(np.inner(k,x)- omega*t[1])) )

    # Plot result if appropriate
    #assert correctResult.all() == result.all(), "Linear2D.height returned %f should be %f" % (result,correctResult)

    try:
        plot(x[0], result, 'b', linewidth=2)
        if showPlots:
            show()
    except:
        pass


def test_Solitary(showPlots=False):
    """ Testing the Linearized 2D interface (phi) propagation. """
    A = amplitude                # amplitude
    k = (2*np.pi/waveLength,0.0,0.0)
    h = L[2]
    omega = np.sqrt(-g[2]*k[0]*np.tanh(k[0]*h))
    period = 2*np.pi/omega

    # Collocation point
    N = 100
    x = [np.linspace(0,L[0],N), 0.0, 0.0]
    t = np.linspace(0,period,N)
    [xx, tt] = np.meshgrid(x,t) 
    
    # Wave Field Object
    waveTest = wm.Solitary(A,omega,k,h,rho_0,rho_1)

    result = waveTest.height(x,t[90])
    #correctResult = np.real( A*np.exp(1j*(k[0]*x[0]- omega*t[1])) )

    # Plot result if appropriate
    #assert correctResult.all() == result.all(), "Linear2D.height returned %f should be %f" % (result,correctResult)

    try:
        plot(x[0], result, 'b', linewidth=2)
        if showPlots:
            show()
    except:
        pass

def test_waveJONSWAP(showPlots=False):
    """ Testing the initialization via JONWAP wave spectrum."""
    # Wave Field Object
    waveTest = wm.waveJONSWAP()

    # Discretization
    Nx = 2**JS.npw1                # number of Fourier modes
    Ny = 2**JS.npw2
    Lx = abs(JS.x1e - JS.x1b)         # domain length
    Ly = abs(JS.y1e - JS.y1b)
    kc_x = 2*np.pi/Lx
    kc_y = 2*np.pi/Ly              # std. wave factors if modes would be integers
    kc_x_modes = Nx*kc_x        
    kc_y_modes = Ny*kc_y        # wave factors for modes obtained from fftfreq()

    k = (2*np.pi/waveLength,0.0,0.0)
    h = L[2]
    omega = np.sqrt(-g[2]*k[0]*np.tanh(k[0]*h))
    period = 2*np.pi/omega

    # Collocation point
    N = 100
    x = [np.linspace(0,L[0],N), np.linspace(0,L[1],N), 0.0]
    t = np.linspace(0,period,N)
    [xx, tt] = np.meshgrid(x,t)

    # Get initial wave field for all x,y at t=0
    result = waveTest.height(x,t[0])

    # Plot result if appropriate
    #assert correctResult.all() == result.all(), "Linear2D.height returned %f should be %f" % (result,correctResult)

    try:
        # Attaching 3D axis to the figure
        fig = plt.figure()
        ax = p3.Axes3D(fig)
        
        x = np.linspace(JS.x1b, JS.x1e, Nx)
        y = np.linspace(JS.y1b,JS.y1e, Ny)
        [xx, yy] = np.meshgrid(x, y)
        
        surf = ax.plot_surface(xx,yy,result,rstride=2,cstride=2, cmap=cm.jet,
                linewidth=0.5, antialiased=False)
        #ax.plot_wireframe(xx,yy,surface, rstride=4, cstride=4)
        plt.show()
    except:
        pass
    
if __name__ == '__main__':
    print "The program name is: ", __name__
    #test_Linear2D(showPlots=True)
    test_WaveGroup(showPlots=True)
    #test_Solitary(showPlots=True)
    #test_waveJONSWAP(showPlots=True)
