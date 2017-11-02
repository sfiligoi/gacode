include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

#DO NOT build superLU by default
default:

all: slu slumt

slu:
	cd SuperLU ; make "OPT=${OPT}";

slumt:
	cd SuperLU_MT ; make "OPT=${OPT}";

clean:
	cd SuperLU ; make clean
	cd SuperLU_MT ; make clean

