import numpy as np
""" Set of wave modules for initializing and
    driving the wavemaker field, that will be
    fed in and constrasted with field data. """


def linear2DHeight(amplitude, omega, k, t, x):
    """ Gives a linearized solution velocity in x-dir to the potential flow
        model in two dimensions (x,y=0,z) for finite depth and infinity and 
        shllow water limits."""
    A = amplitude
    eta = A*np.exp(1j*k*x - omega*t)
    return real(eta)

def linear2DPressure(amplitude, omega, k, t, x, depth, z, rho):
    """ Gives linearized pressured with P_atm = 0 """
    g = 9.81
    p = rho*g*amplitude * np.cosh(k*(z+h))/np.cosh(k*h) * np.exp(1j*k*x - 1j*omega*t)
    return real(p)

def linear2DVelocity_u(amplitude, omega, k, t, x, depth, z):
    """ Defines a linearized solution to the potential flow
        model in two dimensions (x,y=0,z) for finite depth 
        and infinity and shllow water limits."""
    p_atm = 0.0                       # atmospheric pressure
    g = 9.81                          # gravity
    A = amplitude
    h = depth

    # Finite Depth (0 < kh < infty)
    u = (g*k*A / omega) * np.cosh(k*(z+h))/np.cosh(k*h) * np.exp(1j*k*x - 1j*omega*t)  

    # Deep water (kh >> 1)
    # ... TODO
                 
    # Shallow water (kh << 1)
    # ... TODO

    return real(u)


def linear2DVelocity_v(amplitude,omega,k,t,x,depth):
    return v = 0.0


def linear2DVelocity_w(amplitude, omega, k,t, x, depth, z):
    """ Gives a linearized solution velocity in x-dir to the potential flow
        model in two dimensions (x,y=0,z) for finite depth and infinity and 
        shllow water limits."""
    p_atm = 0.0                       # atmospheric pressure                                                                                                                                  
    g = 9.81                          # gravity
    A = amplitude
    h = depth

    # Finite Depth (0 < kh < infty)                                                                                                                                      
    w = -1j * (g*k*A/omega) * np.sinh(k*(z+h))/np.cosh(k*h) * np.exp(1j*k*x - 1j*omega*t)

    # Deep Water (kh >> 1)
    # ... TODO

    # Shallow Water
    # ... TODO

    return real(w)
