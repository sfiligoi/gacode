include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

FC += ${FOMP}

.SUFFIXES : .o .f90

EXTRA_LIBS = \
	${GACODE_ROOT}/gapy/vpro/vpro_lib.a \
	${GACODE_ROOT}/gapy/src/geo_lib.a 

vpro_test: vpro_test.o
	$(FC) -o vpro_test vpro_test.o ${EXTRA_LIBS} ${LMATH}

.f90.o :
	$(FC) $(FFLAGS) $(FMATH) -c $<

clean:
	rm -f *.o vpro_test input.gacode
