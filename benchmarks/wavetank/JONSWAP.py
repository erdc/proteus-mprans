"""
A script for initializing data files of a Quasi-3D wave simulation via
the JONSWAP spectrum and directional distribution.

GENERAL NOTES
=============

PROGRAM:  JONSWAP.py (Python Script File)

AUTHOR:   Matt Malej (matt.malej@erdc.dren.mil)

PURPOSE:  This Python script is designed to set up and initalize the wave
          field via a JONSWAP spectrum and for both directionally confined and
          and directional varied distribution.
          Will serve as input for simulating linear and eventually weakly nonlinear
          random surface water waves in realistic settings. 

          The module imports JONSWAP_p.py (p = physics), where all the input
          parameters are specified. 

INPUT:    Phyical parameters from JONSWAP_p.py

Self NOTE: for FFT your transforming function needs to be periodic on your 
           domain and do NOT make your first and final point in the domain equal!
           ... i.e. for sin(x) do NOT make x=[0,2*pi], but instead x=[dx,2*pi]!

LAST UPDATE: September 10, 2012
"""

def JONSWAP():
    import sys                  # load system and math modules
    from math import *
    import numpy as np          # modules for computing FFT and other numerical stuff
    import JONSAWP_py as JS

    from spectrum import jonswap
    from potential import velocityPotential

    try:
        from matplotlib.pylab import *           # for temporary plots (checking)
    except:
        pass

    # Some Global Variables
    omega_peak = 2.0*pi*JS.fp              # peak angular frequency
    kp = (2.0 * pi * JS.fp)**2 /JS.gv      # peak wavenumber (determinied from disper. rel. with h-->infty)
    neqs = 2                               # number of equations

    # Discretization
    Nx = 2**JS.npw1                # number of Fourier modes
    Ny = 2**JS.npw2
    Lx = abs(JS.x1e - JS.x1b)         # domain length
    Ly = abs(JS.y1e - JS.y1b)
    kc_x = 2*pi/Lx
    kc_y = 2*pi/Ly              # std. wave factors if modes would be integers
    kc_x_modes = Nx*kc_x        
    kc_y_modes = Ny*kc_y        # wave factors for modes obtained from fftfreq()

    modes_x = np.fft.fftfreq(Nx)    # Fourier mode alignments
    modes_y = np.fft.fftfreq(Ny)

    kx = kc_x_modes * modes_x
    ky = kc_y_modes * modes_y

    # Generating mesh grid for 2D spectrum
    [kkx, kky] = np.meshgrid(kx,ky)

    # Instantiate 2D array for the JONSWAP spectrum (complex entries)!
    spec = np.zeros((Nx,Ny), complex)

    # ~ Call jonswap() in spectrum.py module *** Send 1D vectors kx and ky (not kkx nor kky)!
    spectrum = jonswap(spec, JS.Hs, JS.gamma, JS.fp, JS.nspread, JS.thetam, JS.gv, kx, ky)

    # (DISREGARD FOR NOW) - Imposing iniitial filter (same as Trulsen's et al. for mNLS if filter==1)
    if JS.filter == 1:
        spectrum[ 2*floor(kp/kc_x)+1 : (Nx-1)-2*floor(kp/kc_x), 2*floor(kp/kc_x)+1 : (Ny-1)-2*floor(kp/kc_x) ] = 0.0 + 0.0j 
    else:
        pass

    # *** Compute the surface elvation via iFFT2D of the spectrum  --> take REAL part
    #     just in case we have some left over (nonzero) imaginary part of 'spectrum'
    surface = np.real(np.fft.ifft2(spectrum) )

    # *** Compute the velocity potential from linear theory
    velocityPotential = velocity_potential(spectrum, JS.gv, kx, ky)

    # *** Putting the surface(s) into a binary file 
    #
    #    NOTE: Numpy Array has C (row-based) alignemnt and hence needs to be
    #          transposed before writing to a file if these files will be read
    #          in via Fortran90 and its arrays are column-based!  
    # ---> surface = surface.transpose()
    # ---> velocityPotential = velocityPotential.transpose()
    #
    #fout = open("in2.dat", 'wb')                         # binary write mode
    #fwrite(fout, surface.size, surface)                  # surface is REAL 
    #fout.close()
    #
    #fout = open("in2.dat", 'a+b')                               # binary append write mode
    #fwrite(fout, velocityPotential.size, velocityPotential)     # velocityPotential is REAL
    #fout.close()

    try:
        ion()
        figure(1)
        clf()
        plot(kx[0:Nx/2-1,], abs(spectrum[0:Nx/2-1,0]), 'b', linewidth=2)
    except:
        pass

    # Returning surface elvation along the line of wavemakers [x,y]
    location = 10
    return surface[10,:]
