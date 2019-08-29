include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

OBJ=gacode

SRC= 	../vis/gapy_realfluct.f90 \
	../vis/gapy_wheel1.f90 \
	../vis/gapy_wheel2.f90 \
	../vis/gapy_torcut.f90 \
	../vis/gapy_surfpar.f90 \
	../geo/geo.f90 \
	../expro/expro.f90 \
	../expro/expro_util.f90 \
	../expro/expro_pycomm.f90

# Remove possible compilers
F2PYM:=$(filter-out --fcompiler=pg --fcompiler=intelem,$(F2PY))

ifndef F2PY
   all:
else
   all: objs
endif

objs:
	${F2PYM} --quiet --overwrite-signature -m ${OBJ} -h sgn.pyf ${SRC}
	${F2PY}  --quiet -lgomp --f90flags='${FOMP} ${FOPT} ${FMATH}' -c sgn.pyf ${SRC}

	cp *.so ${GACODE_ROOT}/python/

clean:
	rm -f sgn.pyf *~ *.so ${GACODE_ROOT}/python/*.so

exit:
	echo "No f2py"
