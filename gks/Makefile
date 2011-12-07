include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

export EXTRA_LIBS = ${GACODE_ROOT}/tglf/src/tglf_lib.a  \
                    ${GACODE_ROOT}/dataread/src/dataread_lib.a \
                    ${GACODE_ROOT}/shared/math/math_lib.a  \
                    ${LMATH}  \
                    ${NETCDF}  

all: lib
	cd src ; make
lib:
	cd ${GACODE_ROOT}/tglf ; make
	cd ${GACODE_ROOT}/dataread ; make 
	cd ${GACODE_ROOT}/shared/math ; make 

clean:
	cd src ; make clean 

