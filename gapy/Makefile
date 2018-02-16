#---------------------------
# Toplevel makefile for gapy
#---------------------------

include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

export EXTRA_LIBS = \
       ${GACODE_ROOT}/shared/math/math_lib.a 


all:
	cd ${GACODE_ROOT}/shared ; make
	cd src ; make

clean:
	cd src ; make clean

.IGNORE:
