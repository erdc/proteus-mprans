Wavetank (PROTEUS-MPRANS) Sphinx Documentation
==============================================

This is the root directory of the documentation. The documentation is
built from reStructuredText and the python docstrings using sphinx. 

The location of the Makefile should be in your ``doc`` directory under the
root Proteus-MPRANS directory.

To build the html type ``make html``, or to build latex docs ``make latex``,
 or see the Makefile for additional formats. 
The documentation starts in index.rst and is organized as a tree 
using the toctree entry in each .rst file. 

The api, scripts, and test docs are/will be generated from the docstrings 
using the sphinx autodoc extension. See the files in api and tests for examples. 

Note that for autodoc to work the python interpreter must have access to
the files. The scripts and test/problemDescriptions paths are thus
appended in conf.py so that these modules, are accessible to python 
during the sphinx build process.

Additional Notes
================

Remember to export PROTUES_MPRANS ariable, giving your full directory path to
your Proteus-MPRANS installation directory.

For me it is for example (using Bash)::

  export PROTEUS_MPRANS=$HOME/proteus-mprans  #or wherever you installed it

Then::

  cd $PROTEUS_MPRANS/doc

and you are ready to fly!

Also, if you modify the package modules and would like to update the api directory, 
you can build the sphinx api documentation by first invoking.

::

    sphinx-apidoc -o api $PROTUES_MPRANS/doc

Note that you will have to ``rm -rf api`` as Sphinx will not overwrite files by
default. Alternatively, you can add ``-f`` flag to sphinx-apidoc above to force
overwriting of files.

Then you can build html by::

  sphinx-build -b html $PROTEUS_MPRANS/doc _build/

or building LaTex file (where then you can typeset with Latex interpreter) by

::

  sphinx-build -b latex $PROTEUS_MPRANS/doc _build/


Developer Information
======================

The source code, wiki, and issue tracker are on github at
https://github.come/erdc-cm/proteus-mprans. The repository is private and 
the developers' mailing list is http://groups.google.com/group/proteus-dev. 
Both require approval at this time.