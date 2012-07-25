include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

export EXTRA_LIBS = \
	${GACODE_ROOT}/shared/math/math_lib.a \
	${GACODE_ROOT}/shared/legacyread/src/legacyread_lib.a \
	${LMATH} \
	${NETCDF}

all:
	cd ${GACODE_ROOT}/shared/math ; make
	cd ${GACODE_ROOT}/shared/GEO ; make
	cd ${GACODE_ROOT}/shared/EXPRO ; make
	cd ${GACODE_ROOT}/shared/legacyread ; make
	cd locpargen ; make
	cd src ; make

clean:
	cd locpargen ; make clean
	cd src ; make clean
