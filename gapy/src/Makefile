include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

all:
	f2py -m gapy -h sgn.pyf gapy.f90
	f2py -c sgn.pyf gapy.f90 ${GACODE_ROOT}/shared/math/math_lib.a ${LMATH}

clean:
	rm -f gapy.so sgn.pyf *~
