#---------------------------
# Toplevel makefile for gapy
#---------------------------

include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

all:
	cd src  ; make
	cd f2py ; make

clean:
	cd src ; make clean
	cd f2py ; make clean

.IGNORE:
