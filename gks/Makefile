include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

export EXTRA_LIBS = ${GACODE_ROOT}/tglf/src/tglf_lib.a  \
                    ${GACODE_ROOT}/gks/legacyread/legacyread_lib.a \
                    ${GACODE_ROOT}/shared/math/math_lib.a  \
                    ${LMATH}  \
                    ${NETCDF}  

all:
	cd ${GACODE_ROOT}/shared/math ; make
	cd ${GACODE_ROOT}/tglf ; make
	cd legacyread ; make
	cd src ; make

clean:
	cd src ; make clean
