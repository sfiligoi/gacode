include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

OBJ=vis

SRC=vis.f90

# Remove possible compilers
#F2PYM:=$(filter-out --fcompiler=pg --fcompiler=intelem,$(F2PY))

ifndef F2PY
   all:
else
   all: objs
endif

objs:
	${F2PY} --quiet --overwrite-signature -m ${OBJ} -h sgn.pyf ${SRC}
	${F2PY} --quiet -lgomp --f90flags='${FOMP} ${FOPT} ${FMATH}' -c sgn.pyf ${SRC}

clean:
	rm -f py*.so sgn.pyf *~

exit:
	echo "No f2py"
