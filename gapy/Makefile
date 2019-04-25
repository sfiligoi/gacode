#---------------------------
# Toplevel makefile for gapy
#---------------------------

all:
	cd geo  ; make
	cd expro ; make
#	cd f2py ; make
	cd test ; make

pip:
	python setup.py build_ext
	python setup.py sdist

clean:
	cd geo ; make clean
	cd expro ; make clean
	cd f2py ; make clean
	cd test ; make clean
	rm -rf dist build MANIFEST 

.IGNORE:
