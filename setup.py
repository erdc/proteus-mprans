import os
import sys


from distutils.core import setup, Extension
from Cython.Distutils import build_ext
import numpy

## \file setup.py setup.py
#  \brief The python script for building proteus
#
#  Set the DISTUTILS_DEBUG environment variable to print detailed information while setup.py is running.
#

# Get Proteus configuration information

from proteus.config import *

###to turn on debugging in c++
##\todo Finishing cleaning up setup.py/setup.cfg, config.py...
from distutils import sysconfig
cv = sysconfig.get_config_vars()

# Get include flags for building Proteus Extensions

import proteus

try:
    import proteus.util
    proteus_include_dir = proteus.util.get_include_dir()
    proteus_model_kernel = [proteus_include_dir + '/ModelFactory.h', proteus_include_dir + '/CompKernel.h']
except:
    import sys
    sys.stderr.write("Unable to import `get_include_dir` from Proteus, is your Proteus up to date?")
    raise

# Ensure that Proteus namespace init.py has been symbolically linked in
if not os.path.exists(os.path.join('.','proteus','__init__.py')) or not os.path.exists(os.path.join('.','proteus','__init__.pyc')):
    # Symbolically link in Proteus __init__.py file
    proteus_init = proteus.__file__
    # if pyc, use py instead if available
    if proteus_init.endswith('pyc') and os.path.exists(proteus_init[:-1]):
        proteus_init = proteus_init[:-1]
    proteus_init_name = os.path.split(proteus_init)[-1]
    os.symlink(proteus_init, os.path.join('.', 'proteus', proteus_init_name))


setup(name='proteus_mprans',
      version='0.9.0',
      description='Proteus modules for simulating free surface fluid/structure interaction',
      author='Chris Kees',
      author_email='christopher.e.kees@usace.army.mil',
      url='https://proteus.usace.army.mil',
      packages=['proteus.mprans'],
      ext_package='proteus.mprans',
      cmdclass = {'build_ext':build_ext},
      ext_modules=[Extension("cNCLS",["proteus/mprans/cNCLS.pyx"],depends=["proteus/mprans/NCLS.h"], language="c++",
                             include_dirs=[numpy.get_include(), proteus_include_dir]),
                   Extension("cMCorr",["proteus/mprans/cMCorr.pyx"],depends=["proteus/mprans/MCorr.h"], define_macros=[('PROTEUS_LAPACK_H',PROTEUS_LAPACK_H),
                                            ('PROTEUS_LAPACK_INTEGER',PROTEUS_LAPACK_INTEGER),
                                            ('PROTEUS_BLAS_H',PROTEUS_BLAS_H)],language="c++",
                             include_dirs=[numpy.get_include(),proteus_include_dir],
                             library_dirs=[PROTEUS_LAPACK_LIB_DIR,
                                           PROTEUS_BLAS_LIB_DIR],
                             libraries=['m',PROTEUS_LAPACK_LIB,PROTEUS_BLAS_LIB],
                             extra_compile_args=PROTEUS_EXTRA_COMPILE_ARGS,
                             extra_link_args=PROTEUS_EXTRA_LINK_ARGS),
                   Extension("cRANS2P",["proteus/mprans/cRANS2P.pyx"], depends=["proteus/mprans/RANS2P.h"] + proteus_model_kernel,
                             language="c++", include_dirs=[numpy.get_include(), proteus_include_dir]),
                   Extension("cRANS2P2D",["proteus/mprans/cRANS2P2D.pyx"],
                             depends=["proteus/mprans/RANS2P2D.h"] + proteus_model_kernel,
                             language="c++",
                             include_dirs=[numpy.get_include(), proteus_include_dir]),
                   Extension("cRDLS",["proteus/mprans/cRDLS.pyx"],
                             depends=["proteus/mprans/RDLS.h"] + proteus_model_kernel,
                             language="c++",
                             include_dirs=[numpy.get_include(), proteus_include_dir]),
                   Extension("cVOF",["proteus/mprans/cVOF.pyx"],
                             depends=["proteus/mprans/VOF.h"] + proteus_model_kernel,
                             language="c++",
                             include_dirs=[numpy.get_include(), proteus_include_dir]),
                   Extension("cMoveMesh",["proteus/mprans/cMoveMesh.pyx"],
                             depends=["proteus/mprans/MoveMesh.h"] + proteus_model_kernel,
                             language="c++",
                             include_dirs=[numpy.get_include(), proteus_include_dir]),
                   Extension("cSW2D",["proteus/mprans/cSW2D.pyx"],
                             depends=["proteus/mprans/SW2D.h"] + proteus_model_kernel,
                             language="c++",
                             include_dirs=[numpy.get_include(), proteus_include_dir],
                             extra_compile_args=PROTEUS_EXTRA_COMPILE_ARGS+['-g'],
                             extra_link_args=PROTEUS_EXTRA_LINK_ARGS+['-g']),
                   Extension("cSW2DCV",["proteus/mprans/cSW2DCV.pyx"],
                             depends=["proteus/mprans/SW2DCV.h"] + proteus_model_kernel,
                             language="c++",
                             include_dirs=[numpy.get_include(), proteus_include_dir],
                             extra_compile_args=PROTEUS_EXTRA_COMPILE_ARGS+['-g'],
                             extra_link_args=PROTEUS_EXTRA_LINK_ARGS+['-g']),
                   Extension("cKappa",["proteus/mprans/cKappa.pyx"],
                             depends=["proteus/mprans/Kappa.h"] + proteus_model_kernel,
                             language="c++",
                             include_dirs=[numpy.get_include(), proteus_include_dir]),
                   Extension("cKappa2D",["proteus/mprans/cKappa2D.pyx"],
                             depends=["proteus/mprans/Kappa2D.h"] + proteus_model_kernel,
                             language="c++",
                             include_dirs=[numpy.get_include(), proteus_include_dir]),
                   Extension("cDissipation",["proteus/mprans/cDissipation.pyx"],
                             depends=["proteus/mprans/Dissipation.h"] + proteus_model_kernel,
                             language="c++",
                             include_dirs=[numpy.get_include(), proteus_include_dir]),
                   Extension("cDissipation2D",["proteus/mprans/cDissipation2D.pyx"],
                             depends=["proteus/mprans/Dissipation2D.h"] + proteus_model_kernel,
                             language="c++",
                             include_dirs=[numpy.get_include(), proteus_include_dir]),
                   ],
      requires=['numpy']
      )
