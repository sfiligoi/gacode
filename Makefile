#---------------------------
# Toplevel makefile for gkcoll
#---------------------------

include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

export EXTRA_LIBS = \
       ${GACODE_ROOT}/shared/EXPRO/EXPRO_lib.a \
       ${GACODE_ROOT}/shared/math/math_lib.a \
       ${GACODE_ROOT}/shared/GEO/GEO_lib.a \
       ${GACODE_ROOT}/shared/UMFPACK/UMFPACK_lib.a \
       ${LMATH}

all:
	cd ${GACODE_ROOT}/shared ; make
	cd src ; make

clean:
	cd ${GACODE}/modules ; rm -f gkcoll*.mod
	cd src ; make clean

.IGNORE:
