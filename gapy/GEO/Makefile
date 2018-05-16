include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

LLIB=GEO_lib

OBJECTS = geo.o

.SUFFIXES : .o .f90

all: geo_lib.a test.o
	$(FC) -o test test.o geo_lib.a

geo_lib.a: $(OBJECTS)
	$(ARCH) geo_lib.a $(OBJECTS)

.f90.o :
	$(FC) $(FFLAGS) $(FMATH) -c $<

clean:
	rm -f *.o *.a *.mod test geo.so sgn.pyf

py:
	${F2PY} --quiet -m geo -h sgn.pyf geo.f90
	${F2PY} --quiet -lgomp --f90flags='${FOMP} ${FOPT}' -c sgn.pyf geo.f90
