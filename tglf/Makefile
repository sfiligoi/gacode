include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

harvest_lib = ${GACODE_ROOT}/shared/harvest_client/libharvest.a

export EXTRA_LIBS = \
	${harvest_lib} \
	${GACODE_ROOT}/shared/math/math_lib.a

all: ${harvest_lib}
	cd src ; make

clean:
	cd src ; make clean 

${harvest_lib}:
	cd ${GACODE_ROOT}/shared/harvest_client ; make
