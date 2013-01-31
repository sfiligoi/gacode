include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

export EXTRA_LIBS = \
	${GACODE_ROOT}/shared/math/math_lib.a \
	${GACODE_ROOT}/shared/EXPRO/EXPRO_lib.a \
	${GACODE_ROOT}/shared/GEO/GEO_lib.a \
	${LMATH} \
	${NETCDF}

all:
	cd ${GACODE_ROOT}/shared/math ; make
	cd ${GACODE_ROOT}/shared/GEO ; make
	cd ${GACODE_ROOT}/shared/EXPRO ; make
	cd locpargen ; make
	cd src ; make
	cd exprotools ; make

clean:
	cd locpargen ; make clean
	cd src ; make clean
	cd exprotools ; make clean
