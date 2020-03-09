#--------------------------------------
# Toplevel makefile for f2py components
#---------------------------------------

all:
	cd geo ; make
	cd expro ; make

pip:
	python ./setup.py build
	python ./setup.py sdist

so:
	python ./setup.py install --home=${GACODE_ROOT}/python
	mv ${GACODE_ROOT}/python/lib/python/* ${GACODE_ROOT}/python/

clean:
	cd geo ; make clean
	cd expro ; make clean
	cd module ; make clean
	cd test ; make clean
	rm -rf dist build MANIFEST sgn.pyf
	cd ${GACODE_ROOT}/python ; rm -rf gacode *-info *.so

.IGNORE:

