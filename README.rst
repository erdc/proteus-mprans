Proteus-MPRANS
======================================================

Proteus (http://proteus.usace.army.mil) is a Python package for
rapidly developing computer models and numerical methods.

Installation
=============

First install the base proteusModule and set the related enviroment setttings.
Then set the following environment variables and type 'make' in the root directory.

For easy pre/postprocessing and general convenience, set the following enviroment variables. 

:: 

  setenv PROTEUS_MPRANS $HOME/proteus_mprans   #or wherever you put it
  setenv PATH $PATH:$PROTEUS_MPRANS/scripts

Benchmarks
===========

The directory benchmarks contains different directories for each benchmark.
Consult the README in each subdirectory to get relevant information on how 
to preprocess, run and postprocess each benchmark case. 
  
Developer Information
======================

The source code, wiki, and issue tracker are on github at
https://github.come/erdc/proteus. The developers' mailing list is
http://groups.google.com/group/proteus-dev. Both require approval at
this time.
