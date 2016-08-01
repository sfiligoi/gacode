include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

harvest_lib = ${GACODE_ROOT}/shared/harvest_client/libharvest.a
ifdef FANN_ROOT
brainfuse_lib = ${GACODE_ROOT}/shared/neural/libbrainfuse.a
else
brainfuse_lib =
endif

export EXTRA_LIBS = \
	${harvest_lib} \
	${brainfuse_lib} \
	${GACODE_ROOT}/shared/math/math_lib.a

all:
	cd ../shared ; make
	cd src ; make

clean:
	cd src ; make clean 
