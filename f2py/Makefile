#--------------------------------------
# Toplevel makefile for f2py components
#---------------------------------------

all:
	cd geo ; make
	cd expro ; make
#	make so

pip:
	cp -r ../python/cgyro .
	python ./setup.py build
	python ./setup.py sdist

so:
	cd module ; make

clean:
	cd geo ; make clean
	cd expro ; make clean
	cd module ; make clean
	cd test ; make clean
	rm -rf dist build MANIFEST sgn.pyf
	rm -rf cgyro

.IGNORE:

