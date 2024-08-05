#----------------------------
# Toplevel makefile for tgyro
#----------------------------

include ../platform/build/make.inc.${GACODE_PLATFORM}
include Makefile.extra
       ${GACODE_ROOT}/shared/landau/landau.o \
       ${GACODE_ROOT}/shared/landau/gyrotransformation.o

all:
	gacode_getversion > .VERSION
	cd ${GACODE_ROOT}/shared ; make
	cd ${GACODE_ROOT}/f2py ; make
	cd ${GACODE_ROOT}/neo ; make
	cd ${GACODE_ROOT}/tglf ; make
	cd ${GACODE_ROOT}/cgyro ; make
	cd src ; make

clean:
	cd src ; make clean
	rm -rf .VERSION

.IGNORE:
