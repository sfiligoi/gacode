include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

DESTLIB=${GACODE_ROOT}/python/gapy.so

all:
	f2py -m gapy -h sgn.pyf gapy_math.f90 gapy_fluct.f90
	f2py -lgomp --f90flags='${FOMP}' -c sgn.pyf gapy_math.f90 gapy_fluct.f90 ${GACODE_ROOT}/shared/math/math_lib.a ${LMATH}
	ln -s ${GACODE_ROOT}/gapy/src/gapy.so ${DESTLIB}

clean:
	rm -f gapy.so ${DESTLIB} sgn.pyf *~
