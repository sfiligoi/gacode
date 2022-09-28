include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

FC += ${FOMP} ${FACC}

.SUFFIXES : .o .f90 .F90

all: cgyro_test_alltoall

cgyro_test_alltoall: cgyro_test_alltoall.F90
	$(FC) $(FMATH) $(FFLAGS) -o cgyro_test_alltoall cgyro_test_alltoall.F90

clean:
	rm -f cgyro_test_alltoall *.o

