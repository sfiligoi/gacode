include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

FC += ${FOMP} ${FACC}

.SUFFIXES : .o .f90 .F90

OBJECTS = cgyro_test_alltoall cgyro_test_alltoall_nosplit

all: $(OBJECTS)

cgyro_test_alltoall: cgyro_test_alltoall.F90
	$(FC) $(FMATH) $(FFLAGS) -D ALLTOALLBUFMUL=1024 -o $@ $<

cgyro_test_alltoall_nosplit: cgyro_test_alltoall.F90
	$(FC) $(FMATH) $(FFLAGS) -DDISABLE_MPI_SPLIT=1 -D ALLTOALLBUFMUL=1024 -o $@ $<

clean:
	rm -f $(OBJECTS) *.o

