#--------------------------------------
# Toplevel makefile for f2py components
#---------------------------------------

VERSION := $(shell cat pygacode/version)

all:
	cd geo ; make
	cd expro ; make

pip:
	python ./setup.py build
	python ./setup.py sdist

vis:
	cd vis ; python -m numpy.f2py -c -m vis vis.f90 --f90flags="-fopenmp" -lgomp

pip-release:
	python -m twine upload --repository pypi dist/pygacode-$(VERSION).tar.gz

conda:
	python ./setup.py conda

clean:
	cd geo ; make clean
	cd expro ; make clean
	cd module ; make clean
	cd test ; make clean
	cd vis ; rm *.so

.IGNORE:

