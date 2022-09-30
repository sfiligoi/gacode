include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

FC += ${FOMP} ${FACC}

.SUFFIXES : .o .f90 .F90

OBJECTS = cgyro_test_alltoall cgyro_test_alltoall_nosplit \
	  cgyro_test_allreduce_4 cgyro_test_allreduce_6 cgyro_test_allreduce_32

all: $(OBJECTS)

cgyro_test_alltoall: cgyro_test_alltoall.F90
	$(FC) $(FMATH) $(FFLAGS) -D ALLTOALLBUFMUL=1024 -o $@ $<

cgyro_test_alltoall_nosplit: cgyro_test_alltoall.F90
	$(FC) $(FMATH) $(FFLAGS) -DDISABLE_MPI_SPLIT=1 -D ALLTOALLBUFMUL=1024 -o $@ $<

cgyro_test_allreduce_4: cgyro_test_allreduce.F90
	$(FC) $(FMATH) $(FFLAGS) -D ALLREDUCEBUFMUL=256 -D ALLTOALLBUFMUL=64 -D ALLTOALLNUM=4 -o $@ $<

cgyro_test_allreduce_6: cgyro_test_allreduce.F90
	$(FC) $(FMATH) $(FFLAGS) -D ALLREDUCEBUFMUL=256 -D ALLTOALLBUFMUL=64 -D ALLTOALLNUM=6 -o $@ $<

cgyro_test_allreduce_32: cgyro_test_allreduce.F90
	$(FC) $(FMATH) $(FFLAGS) -D ALLREDUCEBUFMUL=256 -D ALLTOALLBUFMUL=64 -D ALLTOALLNUM=32 -o $@ $<

clean:
	rm -f $(OBJECTS) *.o

