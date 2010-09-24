from distutils.core import setup, Extension
import numpy

## \file setup.py setup.py
#  \brief The python script for building proteus
#
#  Set the DISTUTILS_DEBUG environment variable to print detailed information while setup.py is running.
#
try:
    from config import *
except:
    raise RuntimeError("Missing or invalid config.py file. See proteusConfig for examples")

###to turn on debugging in c++
##\todo Finishing cleaning up setup.py/setup.cfg, config.py...
from distutils import sysconfig
cv = sysconfig.get_config_vars()
cv["OPT"] = cv["OPT"].replace("-DNDEBUG","-DDEBUG")
cv["OPT"] = cv["OPT"].replace("-O3","-g")
cv["CFLAGS"] = cv["CFLAGS"].replace("-DNDEBUG","-DDEBUG")
cv["CFLAGS"] = cv["CFLAGS"].replace("-O3","-g")

setup(name='proteus_fs3d',
      version='0.9.0',
      description='Python tools for two-phase incompressible flows',
      author='Chris Kees',
      author_email='christopher.e.kees@usace.army.mil',
      url='https://adh.usace.army.mil/proteus',
      package_dir={'proteus.fs3d':'.'},
      packages=['proteus.fs3d'],
      py_modules=['proteus.fs3d.RANS2PV2',
                  'proteus.fs3d.NCLSV2',
                  'proteus.fs3d.VOFV2',
                  'proteus.fs3d.RDLSV2',
                  'proteus.fs3d.MCORRV2'],
      ext_modules=[Extension('proteus.fs3d.cRANS2PV2',
                             ['cRANS2PV2Module.cpp','RANS2PV2.cpp'],
                             define_macros=[('PROTEUS_SUPERLU_H',PROTEUS_SUPERLU_H)],
                             include_dirs=[numpy.get_include(),PROTEUS_INCLUDE_DIR,os.getenv('PROTEUS')+'/proteusModule/include',os.getenv('PROTEUS')+'/proteusModule/proteus',
                                           PROTEUS_SUPERLU_INCLUDE_DIR],
                             libraries=['m'],
                             extra_compile_args=PROTEUS_EXTRA_COMPILE_ARGS,
                             extra_link_args=PROTEUS_EXTRA_LINK_ARGS),
                   Extension('proteus.fs3d.cNCLSV2',
                             ['cNCLSV2Module.cpp','NCLSV2.cpp'],
                             define_macros=[('PROTEUS_SUPERLU_H',PROTEUS_SUPERLU_H)],
                             include_dirs=[numpy.get_include(),PROTEUS_INCLUDE_DIR,os.getenv('PROTEUS')+'/proteusModule/include',os.getenv('PROTEUS')+'/proteusModule/proteus',
                                           PROTEUS_SUPERLU_INCLUDE_DIR],
                             libraries=['m'],
                             extra_compile_args=PROTEUS_EXTRA_COMPILE_ARGS,
                             extra_link_args=PROTEUS_EXTRA_LINK_ARGS),
                   Extension('proteus.fs3d.cVOFV2',
                             ['cVOFV2Module.cpp','VOFV2.cpp'],
                             define_macros=[('PROTEUS_SUPERLU_H',PROTEUS_SUPERLU_H)],
                             include_dirs=[numpy.get_include(),PROTEUS_INCLUDE_DIR,os.getenv('PROTEUS')+'/proteusModule/include',os.getenv('PROTEUS')+'/proteusModule/proteus',
                                           PROTEUS_SUPERLU_INCLUDE_DIR],
                             libraries=['m'],
                             extra_compile_args=PROTEUS_EXTRA_COMPILE_ARGS,
                             extra_link_args=PROTEUS_EXTRA_LINK_ARGS),
                   Extension('proteus.fs3d.cRDLSV2',
                             ['cRDLSV2Module.cpp','RDLSV2.cpp'],
                             define_macros=[('PROTEUS_SUPERLU_H',PROTEUS_SUPERLU_H)],
                             include_dirs=[numpy.get_include(),PROTEUS_INCLUDE_DIR,os.getenv('PROTEUS')+'/proteusModule/include',os.getenv('PROTEUS')+'/proteusModule/proteus',
                                           PROTEUS_SUPERLU_INCLUDE_DIR],
                             libraries=['m'],
                             extra_compile_args=PROTEUS_EXTRA_COMPILE_ARGS,
                             extra_link_args=PROTEUS_EXTRA_LINK_ARGS),
                   Extension('proteus.fs3d.cMCORRV2',
                             ['cMCORRV2Module.cpp','MCORRV2.cpp'],
                             define_macros=[('PROTEUS_SUPERLU_H',PROTEUS_SUPERLU_H)],
                             include_dirs=[numpy.get_include(),PROTEUS_INCLUDE_DIR,os.getenv('PROTEUS')+'/proteusModule/include',os.getenv('PROTEUS')+'/proteusModule/proteus',
                                           PROTEUS_SUPERLU_INCLUDE_DIR],
                             libraries=['m'],
                             extra_compile_args=PROTEUS_EXTRA_COMPILE_ARGS,
                             extra_link_args=PROTEUS_EXTRA_LINK_ARGS),
                   ],
      requires=['numpy']
      )
