import waveModules_Matt as wm
import JONSWAP_p as JS
import numpy as np

try:
    #from enthought.mayavi import mlab      # for 3D plots 
    import animation                        # in ~/proteus-mprans/.../wavetank/ as my version of matplotlib doesn't have it
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import mpl_toolkits.mplot3d.axes3d as p3
except:
    pass

""" Testing different analytic solution for wavetank benchmark case."""

g = (0.0,0.0,-9.81)        # gravity
L = (20.0,0.25,1.0)         # L[0]=20.0 ...tank dimensions

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

A = amplitude                # amplitude
k = (2*np.pi/waveLength,0.0,0.0)
h = 0.5*L[2]
omega = np.sqrt(-g[2]*k[0]*np.tanh(k[0]*h))
period = 2*np.pi/omega

#see nose docs for more complex testing

def test_Linear2D(showPlots=False):
    """ Testing the Linearized 2D interface (phi) propagation. """

    # Collocation point
    N = 200
    x = [np.linspace(0,L[0],N), 0.0, 0.0]
    t = np.linspace(0,10*period,N)
    #[xx, tt] = np.meshgrid(x,t) 
    
    # Wave Field Object
    waveTest = wm.Linear2D(A,omega,k,h,rho_0,rho_1)

    # Plot result if appropriate
    #assert correctResult.all() == result.all(), "Linear2D.height returned %f should be %f" % (result,correctResult)

    fig = plt.figure()
    y = waveTest.height(x,t[0]) 
    line, = plt.plot(x[0],y)
    plt.axis((0.0, L[0], 0, L[2]))

    def animate(i):
        line.set_ydata(waveTest.height(x,t[i]))  # update the data
        return line,

    def init():
        line.set_ydata(np.ma.array(x[0], mask=True))
        return line,
        
    try:
        anim = animation.FuncAnimation(fig, animate, np.arange(1, 200), init_func=init,
            interval=25, blit=True)    

        #anim.save('monochromatic_wave_BC_animation.mov', fps=20)
        plt.show()
    except:
        pass

def test_WaveGroup(showPlots=False):
    """ Testing the Linearized 2D interface (phi) propagation. """

    # Collocation point
    N = 1000
    x = [np.linspace(0,L[0],N), 0.0, 0.0]
    t = np.linspace(0,10*period,N)
    [xx, tt] = np.meshgrid(x,t) 
    
    # Wave Field Object
    waveTest = wm.WaveGroup(A,omega,k,h,rho_0,rho_1)

    # Plot result if appropriate
    #assert correctResult.all() == result.all(), "Linear2D.height returned %f should be %f" % (result,correctResult)

    fig = plt.figure()
    y = waveTest.height(x,t[0])
    line, = plt.plot(x[0],y)
    plt.axis((0.0, L[0], 0.0, L[2]))

    def animate(i):
        line.set_ydata(waveTest.height(x,t[i]))  # update the data
        return line,

    def init():
        line.set_ydata(np.ma.array(x[0], mask=True))
        return line,
        
    try:
        anim = animation.FuncAnimation(fig, animate, np.arange(1, 1000), init_func=init,
            interval=2, blit=True)    

        #anim.save('modulated_waveGroup_BC_animation.mov', fps=30)
        plt.show()
    except:
        pass


def test_Solitary(showPlots=False):
    """ Testing the Linearized 2D interface (phi) propagation. """

    # Collocation point
    N = 500
    x = [np.linspace(0,L[0],N), 0.0, 0.0]
    t = np.linspace(0,10*period,N)
    [xx, tt] = np.meshgrid(x,t) 
    
    # Wave Field Object
    waveTest = wm.Solitary(2.0*A,omega,k,h,rho_0,rho_1)

    fig = plt.figure()
    y = waveTest.height(x,t[0])
    line, = plt.plot(x[0],y)
    plt.axis((0.0, L[0], 0, L[2]))

    def animate(i):
        line.set_ydata(waveTest.height(x,t[i]))  # update the data
        return line,

    def init():
        line.set_ydata(np.ma.array(x[0], mask=True))
        return line,
        
    try:
        anim = animation.FuncAnimation(fig, animate, np.arange(1, 500), init_func=init,
            interval=20, blit=True)    

        #anim.save('solitary_wave_BC_animation.mov', fps=30)
        plt.show()
    except:
        pass


def test_StokesWave(showPlots=False):
    """ Testing the 2nd order nonlinear Stokes Wave solution. """

    # Collocation point
    N = 500
    x = [np.linspace(0,L[0],N), 0.0, 0.0]
    t = np.linspace(0,10*period,N)
    [xx, tt] = np.meshgrid(x,t) 
    
    # Wave Field Object
    waveTest = wm.StokesWave(A,omega,k,h,rho_0,rho_1)

    fig = plt.figure()
    y = waveTest.height(x,t[0])
    line, = plt.plot(x[0],y)
    plt.axis((0.0, L[0], 0, L[2]))

    def animate(i):
        line.set_ydata(waveTest.height(x,t[i]))  # update the data
        return line,

    def init():
        line.set_ydata(np.ma.array(x[0], mask=True))
        return line,
        
    try:
        anim = animation.FuncAnimation(fig, animate, np.arange(1, 500), init_func=init,
            interval=25, blit=True)    
        #anim.save('stokes_wave_BC_animation.mov', fps=30)
        plt.show()
    except:
        pass


def test_waveJONSWAP(showPlots=False):
    """ Testing the initialization via JONWAP wave spectrum."""
    # Wave Field Object
    waveTest = wm.waveJONSWAP(A,h,L)

    # Discretization
    Nx = 2**JS.npw1                # number of Fourier modes
    Ny = 2**JS.npw2
    kc_x = 2*np.pi/L[0]
    kc_y = 2*np.pi/L[1]              # std. wave factors if modes would be integers
    kc_x_modes = Nx*kc_x        
    kc_y_modes = Ny*kc_y        # wave factors for modes obtained from fftfreq()

    # Collocation point
    N = 100
    x = [np.linspace(0,L[0],Nx), np.linspace(0,L[1],Ny), 0.0]
    t = np.linspace(0,10*period,N)
    [xx, yy] = np.meshgrid(x[0],x[1])

    # Get initial wave field for all x,y at t=0
    result = np.transpose(waveTest.height(x,t[0]))

    try:
        # Attaching 3D axis to the figure
        fig = plt.figure()
        ax = p3.Axes3D(fig)
        
        surf = ax.plot_surface(xx,yy,result,rstride=2,cstride=2, cmap=cm.jet,
                linewidth=0.5, antialiased=False)
        #ax.plot_wireframe(xx,yy,surface, rstride=4, cstride=4)
        plt.show()
    except:
        pass
    
if __name__ == '__main__':
    print "The program name is: ", __name__
    test_Linear2D(showPlots=True)
    test_WaveGroup(showPlots=True)
    test_Solitary(showPlots=True)
    test_StokesWave(showPlots=True)
    test_waveJONSWAP(showPlots=True)
