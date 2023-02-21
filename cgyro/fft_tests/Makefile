include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

ifdef FFTW_INC
   FFLAGS += -I${FFTW_INC}
endif

FC += ${FOMP} ${FACC}

.SUFFIXES : .o .f90 .F90

all: cgyro_test_fft cgyro_test_fft_mpi

cgyro_test_fft: cgyro_test_fft.o $(EXTRA_LIBS)
	$(FC) $(FFLAGS) -o cgyro_test_fft cgyro_test_fft.o $(EXTRA_LIBS) $(LMATH)

cgyro_test_fft_mpi: cgyro_test_fft_mpi.o $(EXTRA_LIBS)
	$(FC) $(FFLAGS) -o cgyro_test_fft_mpi cgyro_test_fft_mpi.o $(EXTRA_LIBS) $(LMATH)

cgyro_test_fft.o: cgyro_test_fft.F90
	$(FC) $(FMATH) $(FFLAGS) -o cgyro_test_fft.o -c cgyro_test_fft.F90

cgyro_test_fft_mpi.o: cgyro_test_fft.F90
	$(FC) $(FMATH) $(FFLAGS) -DTEST_MPI -o cgyro_test_fft_mpi.o -c cgyro_test_fft.F90

clean:
	rm -rf *.o cgyro_test_fft cgyro_test_fft_mpi
