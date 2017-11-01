include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

all:
	cd SuperLU ; make "OPT=${OPT}"
	cd SuperLU_MT ; make "OPT=${OPT}"

clean:
	cd SuperLU ; make clean
	cd SuperLU_MT ; make clean

