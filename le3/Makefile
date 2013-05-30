#---------------------------
# Toplevel makefile for 3d
#---------------------------

include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

export EXTRA_LIBS = \
       ${GACODE_ROOT}/le3/minpack/minpack_lib.a \
       ${LMATH}

all:
	gacode_getversion > .VERSION
	cd minpack ; make
	cd src ; make

clean:
	cd minpack ; make clean 
	cd src ; make clean
	rm -f .VERSION
.IGNORE:
