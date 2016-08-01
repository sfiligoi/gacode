include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

export EXTRA_LIBS = ${GACODE_ROOT}/shared/math/math_lib.a

all:
	cd ../shared ; make
	cd src ; make

clean:
	cd src ; make clean 
