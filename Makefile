default: install_proteus install 

PROTEUS ?= $(shellpwd)
PROTEUS_ARCH ?= $(shell python -c "import sys; print sys.platform")
PROTEUS_PREFIX ?= ${PROTEUS}/${PROTEUS_ARCH}
PROTEUS_PYTHON ?= ${PROTEUS_PREFIX}/bin/python

install_proteus:
	cd ${PROTEUS} && make
	echo ${PWD}

install:
	echo ${PROTEUS_PYTHON}
	echo ${PROTEUS_PREFIX}/bin/python
	${PROTEUS_PYTHON} setup.py install

clean:
	${PROTEUS_PYTHON} setup.py clean
	cd scripts && make clean

cleaner: clean
	touch mprans/*.h

mpkg:
	bdist_mpkg --readme=README --license=COPYING --welcome=INTRODUCTION

source:
	${PROTEUS_PYTHON} setup.py sdist

binary:
	${PROTEUS_PYTHON} setup.py bdist_dumb

rpm:
	${PROTEUS_PYTHON} setup.py bdist_rpm

wininst:
	${PROTEUS_PYTHON} setup.py bdist_wininst
