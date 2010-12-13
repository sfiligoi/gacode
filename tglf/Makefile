#---------------------------
# Toplevel makefile for tglf
#---------------------------

include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

export EXTRA_LIBS = ${LMATH}

all:
	cd src ; make

clean:
	cd ../modules ; rm -f tglf*.mod
	cd src ; make clean 

.IGNORE:
