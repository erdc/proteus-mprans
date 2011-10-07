default: install_proteus install 

install_proteus:
	cd ${PROTEUS} && make

install:
	${PROTEUS_PYTHON} setup.py install
	cd scripts && make

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
