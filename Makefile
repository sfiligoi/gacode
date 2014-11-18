#----------------------------
# Toplevel makefile for cgyro
#----------------------------

include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

export EXTRA_LIBS = \
       ${GACODE_ROOT}/shared/GEO/GEO_lib.a \
       ${GACODE_ROOT}/shared/math/math_lib.a \
       ${LMATH}

all:
	cd ${GACODE_ROOT}/shared/math ; make
	cd ${GACODE_ROOT}/shared/GEO ; make
	cd src ; make

clean:
	cd ${GACODE}/modules ; rm -f cgyro*.mod
	cd src ; make clean

.IGNORE:
