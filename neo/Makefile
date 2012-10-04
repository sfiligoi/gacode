#---------------------------
# Toplevel makefile for neo
#---------------------------

include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

export EXTRA_LIBS = \
       ${GACODE_ROOT}/shared/EXPRO/EXPRO_lib.a \
       ${GACODE_ROOT}/shared/math/math_lib.a \
       ${GACODE_ROOT}/shared/GEO/GEO_lib.a \
       ${GACODE_ROOT}/shared/UMFPACK/UMFPACK_lib.a \
       ${GACODE_ROOT}/neo/nclass/nclass_lib.a \
       ${LMATH}

all:
	gacode_getversion > .VERSION
	cd ${GACODE_ROOT}/shared ; make
	cd nclass ; make
	cd src ; make
	cd tools/vgen ; make

clean:
	cd ../modules ; rm -f neo*.mod .VERSION
	cd src ; make clean
	cd tools/vgen ; make clean 
	cd nclass ; make clean 

.IGNORE:
