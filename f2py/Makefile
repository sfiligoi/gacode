#--------------------------------------
# Toplevel makefile for f2py components
#---------------------------------------

all:
	cd geo ; make
	cd expro ; make

pip:
	python setup.py build_ext
	python setup.py sdist

so:
	cd module ; make

clean:
	cd geo ; make clean
	cd expro ; make clean
	cd module ; make clean
	cd test ; make clean
	rm -rf dist build MANIFEST sgn.pyf

.IGNORE:

