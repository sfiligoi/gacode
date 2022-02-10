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

pip-test-release:
	cd dist; pip install pygacode-$(VERSION).tar.gz

pip-release:
	python -m twine upload --repository pypi dist/pygacode-$(VERSION).tar.gz

conda:
	python ./setup.py conda

so:
	python ./setup.py install --home=${GACODE_ROOT}/f2py/pygacode/
	cd pygacode ; ln -s lib/python/gacode*.so .

clean:
	cd geo ; make clean
	cd expro ; make clean
	cd module ; make clean
	cd test ; make clean
	rm -rf dist build MANIFEST sgn.pyf pygacode/lib *-info
	cd ${GACODE_ROOT}/f2py/pygacode ; rm -rf gacode *-info *.so

.IGNORE:

