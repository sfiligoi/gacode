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

EXEC = cgyro_test_fft

OBJECTS = cgyro_test_fft.o 

.SUFFIXES : .o .f90 .F90

all: $(OBJECTS) $(EXTRA_LIBS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS) $(EXTRA_LIBS) $(LMATH)

.f90.o .F90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -rf *.o $(EXEC)

