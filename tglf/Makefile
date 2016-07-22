include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

harvest_lib = ${GACODE_ROOT}/shared/harvest_client/libharvest.a
brainfuse_lib = ${GACODE_ROOT}/shared/neural/libbrainfuse.a

export EXTRA_LIBS = \
	${harvest_lib} \
	${brainfuse_lib} \
	${GACODE_ROOT}/shared/math/math_lib.a

all: ${harvest_lib} ${brainfuse_lib}
	cd ../shared/neural/; make
	cd src ; make

clean:
	cd src ; make clean 

${harvest_lib}:
	cd ${GACODE_ROOT}/shared/harvest_client ; make

${brainfuse_lib}:
	cd ${GACODE_ROOT}/shared/neural ; make
