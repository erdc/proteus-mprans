import numpy as np
from math import pi
import JONSWAP_p as JS
""" Set of wave modules for initializing and
    driving the wavemaker field, that will be
    fed in and constrasted with wave flume 
    experimental data."""

class Linear2D:
    """
    A class for linearized solutions of 2D flow (1D + surface) for
    travelling progressive waves ~ exp(kx-wt)
    
    .. todo:: Finish the docs

    An equation
    
    .. math:: 
        
        \Phi_t = -g \zeta + ...

    More text, inline math :math:`x^3` 
    """
    #see this url for some useful restructured text directives
    #
    #http://matplotlib.sourceforge.net/devel/documenting_mpl.html#formatting-mpl-docs
    def __init__(self,amplitude,omega,k,depth,rho_0,rho_1):
        self.A = amplitude
        self.omega = omega
        self.k = k
        self.h = depth
        self.rho_0 = rho_0      # density of water
        self.rho_1 = rho_1      # density of air

    def height(self,x,t):
        """ Gives a linearized solution for the air-water interface to the 
            potential flow model in two dimensions (x,y,z=eta) for finite depth."""
        eta = self.A*np.exp(1j*(self.k[0]*x[0] - self.omega*t))
        # ~ NOTE: x[0] is a vector here!
        return np.real(eta)

    def pressure(self,x,t):
        """ Gives linearized pressured with P_atm = 0 """
        g = (0.0,0.0,-9.81)
        p = self.rho_0*(-g[2])*self.A* np.cosh(self.k[0]*(x[2]+self.h))/np.cosh(self.k[0]*self.h) * np.exp(1j*(self.k[0]*x[0] - self.omega*t))
        return np.real(p)

    def velocity_u(self,x,t):
        """ Defines a linearized solution to the potential flow
            model in two dimensions (x,y,z) for finite depth,
            as well as, deep and shllow water limits.

            .. todo:: implement deep & shallow water limits."""
        g = (0.0,0.0,-9.81)                          # gravity

        # Finite Depth (0 < kh < infty)
        u = (-g[2]*self.k[0]*self.A / self.omega) * np.cosh(self.k[0]*(z+self.h))/np.cosh(self.k[0]*self.h) * np.exp(1j*(self.k[0]*x[0] - self.omega*t))  

        # Deep water (kh >> 1)
        # ... TODO
                 
        # Shallow water (kh << 1)
        # ... TODO

        return np.real(u)


    def velocity_v(self,x,t):
        v = 0.0
        return v

    def velocity_w(self,x,t):
        """ Gives a linearized solution velocity in x-dir to the potential flow
            model in two dimensions (x,y,z) for finite depth, as well as, deep
            and shallow water limits.

            .. todo:: implement deep & shallow water limits."""
        g = (0.0,0.0,-9.81)                          # gravity
        
        # Finite Depth (0 < kh < infty)                                                                                                                                      
        w = -1j * (-g[2]*self.k[0]*self.A/self.omega) * np.sinh(self.k[0]*(x[2]+self.h))/np.cosh(self.k[0]*self.h) * np.exp(1j*(self.k[0]*x[0] - self.omega*t))

        # Deep Water (kh >> 1)
        # ... TODO

        # Shallow Water
        # ... TODO

        return np.real(w)


class WaveGroup:
    """ Class that defines a nearly monochromatic
        wave train/group of waves of same amplitude.
        
        .. todo:: Finish the docs. """

    def __init__(self,amplitude,omega,k,depth,rho_0,rho_1):
        self.A = amplitude
        self.omega = omega
        self.k = k
        self.h = depth
        self.rho_0 = rho_0      # density of water                                      
        self.rho_1 = rho_1      # density of air
        self.N = 2              # number of 

    def height(self,x,t):
        theta =  self.k[0]*x[0] - self.omega*t # ~ NOTE: x[0] is a vector here!
        diff = 0.05
        dtheta = diff*theta
        eta = self.A*np.cos(theta)

        for i in range(self.N):
            eta = eta + self.A*np.sin(theta+(i+1)*dtheta) + self.A*np.sin(theta-(i+1)*dtheta)

        return eta


    def velocity_u(self,x,t):
        """ Defines a linearized solution to the potential flow
            model in two dimensions (x,y,z) for finite depth,
            as well as, deep and shllow water limits, for slowly
            varying regular wavetrains.

            .. todo:: implement deep & shallow water limits. """
        g = (0.0,0.0,-9.81)                          # gravity

        # Finite Depth (0 < kh < infty)
        for i in range(self.N):
            diffPos = (1+diff*(i+1))
            diffNeg = (1-diff*(i+1))
            u = u + (-g[2]*diffPos*self.k[0]*self.A / (diffPos*self.omega)) * np.cosh(diffPos*self.k[0]*(z+self.h))/np.cosh(diffPos*self.k[0]*self.h) * np.exp(1j*diffPos*(self.k[0]*x[0] - self.omega*t)) +\
                (-g[2]*diffNeg*self.k[0]*self.A / (diffNeg*self.omega)) * np.cosh(diffNeg*self.k[0]*(z+self.h))/np.cosh(diffNeg*self.k[0]*self.h) * np.exp(1j*diffNeg*(self.k[0]*x[0] - self.omega*t))
        # Deep water (kh >> 1)
        # ... TODO
                 
        # Shallow water (kh << 1)
        # ... TODO

        return np.real(u)
        # NOTE: implemented based on linearized ideal flow

    def velocity_v(self,x,t):
        v = 0.0
        return v
        # NOTE: you can implement based on linearized ideal flow   


    def velocity_w(self,x,t):
        for i in range(self.N):
            diffPos = (1+diff*(i+1))
            diffNeg = (1-diff*(i+1))        
            w = w + -1j * (-g[2]*diffPos*self.k[0]*self.A/(diffPos*self.omega)) * np.sinh(diffPos*self.k[0]*(x[2]+self.h))/np.cosh(diffPos*self.k[0]*self.h) * np.exp(1j*(diffPos*self.k[0]*x[0] - self.omega*t)) + \
                -1j * (-g[2]*diffNeg*self.k[0]*self.A/self.omega) * np.sinh(self.k[0]*(x[2]+self.h))/np.cosh(self.k[0]*self.h) * np.exp(1j*(self.k[0]*x[0] - self.omega*t))
        
        return w

    def pressure(self,x,t):
        for i in range(self.N):
            diffPos = (1+diff*(i+1))
            diffNeg = (1-diff*(i+1))
            p = p + self.rho_0*(-g[2])*self.A* np.cosh(diffPos*self.k[0]*(x[2]+self.h))/np.cosh(diffPos*self.k[0]*self.h) * np.exp(1j*diffPos*(self.k[0]*x[0] - self.omega*t)) + \
                self.rho_0*(-g[2])*self.A* np.cosh(diffNeg*self.k[0]*(x[2]+self.h))/np.cosh(diffNeg*self.k[0]*self.h) * np.exp(1j*diffNeg*(self.k[0]*x[0] - self.omega*t))
        return np.real(p)
        # NOTE: also implemented on ideal flow via linearized Bernoulli eqn.


class Solitary:
    """ Class that defines a solitary wave profile
        of a constant initial amplitude.
        
        .. todo:: Finish the docs. """
    def __init__(self,amplitude,omega,k,depth,rho_0,rho_1):
        self.A = amplitude
        self.omega = omega
        self.k = k
        self.h = depth
        self.rho_0 = rho_0      # density of water
        self.rho_1 = rho_1      # density of air         
        self.sigma = 4.0        # std. dev.

    def height(self,x,t):
        eta = self.A/np.cosh((self.k[0]*x[0] - self.omega*t)**2 / self.sigma**2)
        # ~ NOTE: x[0] is a vector here!
        return eta

    def velocity_u(self,x,t):
        u = 0.0
        return u
        # NOTE: base it on linearized ideal fluid flow

    def velocity_v(self,x,t):
        v = 0.0
        return v
        # NOTE: base it on linearized ideal fluid flow
    
    def velocity_w(self,x,t):
        w = 0.0
        return w
        # NOTE: base it on ideal fluid flow

    def pressure(self,x,t):
        p = 0.0       # P_atm
        return
        # NOTE: define via Bernoulli eqn.



class waveJONSWAP:
    """ Class that defines a wave field based on realistic wave spectrum.
        The JONSWAP wave spectrum with directional distribution and
        random phases.
        
        .. todo:: Finish the docs. """

    def __init__(self):
        # Based on Linearized Dispersion Relations
        self.omegaPeak = 2.0*pi*JS.fp                  # peak angular frequency
        self.kp = (2.0*pi*JS.fp)**2 /JS.gv             # peak wavenumber 

        # Discretization
        self.Nx = 2**JS.npw1                           # number of Fourier modes
        self.Ny = 2**JS.npw2
        self.Lx = abs(JS.x1e - JS.x1b)                 # domain length
        self.Ly = abs(JS.y1e - JS.y1b)
        self.kc_x = 2*pi/self.Lx
        self.kc_y = 2*pi/self.Ly                        
        self.kc_x_modes = self.Nx*self.kc_x        
        self.kc_y_modes = self.Ny*self.kc_y

        self.surface = np.zeros((self.Nx,self.Ny))      # free surface
        self.u = np.zeros((self.Nx,self.Ny))            # u velocity
        self.v = np.zeros((self.Nx,self.Ny))            # v velocity


    def height(self,x,t):
        [self.surface, self.u, self.v] = self.JONSWAP()
        # ~ NOTE: x[0] is a vector here!
        return self.surface

    def velocity_u(self,x,t):
        """ Returns the velocity at the free-surface (z=surface)
            
            NOTE: to extract velocity at any height z, write down the
                power series and take fft() of for every z
                ... little time consuming, might need a bit of Cython here!

            The velocity potential defined at the free-surace :math:`z=\zeta` is
            given by :math:`\Phi(x,y,t) \equiv \phi(x,y,z=\zeta,t)`
    
            .. math:: 
        
                \phi(x,y,z,t) = \Phi + (z-\zeta)W + \sum ... + \sum ...

            where :math:`W` is the vertical velocity defined at the free-surface
            given through a Dirichlet-to-Neumann operator relation.
        """
        return self.u
        # NOTE: base it on linearized ideal fluid flow

    def velocity_v(self,x,t):
        """ Returns the velocity at the free-surface (z=surface)
            
            NOTE: to extract velocity at any height z, write down the
                power series and take fft() of for every z
                ... little time consuming, might need a bit of Cython here!

            The velocity potential defined at the free-surace :math:`z=\zeta` is
            given by :math:`\Phi(x,y,t) \equiv \phi(x,y,z=\zeta,t)`
    
            .. math:: 
        
                \phi(x,y,z,t) = \Phi + (z-\zeta)W + \sum ... + \sum ...

            where :math:`W` is the vertical velocity defined at the free-surface
            given through a Dirichlet-to-Neumann operator relation.
        """        
        return self.v
        # NOTE: base it on linearized ideal fluid flow

    def JONSWAP(self):
        """Sets a wave field according to JONSWAP ocean wave spectrum."""
        from spectrum import jonswap
        from potential import velocityPotential
       
        modes_x = np.fft.fftfreq(self.Nx)    # Fourier mode alignments
        modes_y = np.fft.fftfreq(self.Ny)

        kx = self.kc_x_modes * modes_x
        ky = self.kc_y_modes * modes_y

        # Generating mesh grid for 2D spectrum
        [kxx, kyy] = np.meshgrid(kx,ky)
    
        # ~ Call jonswap() in spectrum.py module *** Send 1D vectors kx and ky (not kxx nor kyy)!
        spectrum = jonswap(self.Nx, self.Ny, JS.Hs, JS.gamma, JS.fp, JS.nspread, JS.thetam, JS.gv, kx, ky)

        # (DISREGARD FOR NOW) - Imposing iniitial circular zero-pad filter
        if JS.filter == 1:
            spectrum[ 2*floor(self.kp/self.kc_x)+1 : (self.Nx-1)-2*floor(self.kp/self.kc_x), \
                2*floor(selfkp/kc_x)+1 : (self.Ny-1)-2*floor(self.kp/self.kc_x) ] = 0.0 + 0.0j 
        else:
            pass

        # Compute the surface elvation via iFFT2D of the spectrum  --> take REAL part
        # just in case we have some left over (nonzero) imaginary part of 'spectrum'
        surface = np.real(np.fft.ifft2(spectrum) )

        # Compute the velocity potential from linear theory and horizontal velocity
        velPotential, velocity_u, velocity_v = velocityPotential(spectrum, JS.gv, kx, ky)

        # Returning surface elvation along the line of wavemakers [x,y]
        loc = 10
        return surface, velocity_u, velocity_v 

