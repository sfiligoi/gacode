include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

TGLFISNAN = tglf_isnan_generic.o

include ../install/make.ext.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC1=tglf_mpi
EXEC2=tglf
EXEC3=tglf_nn

LLIB=tglf_lib

OBJECTS = ${TGLFISNAN}      \
          tglf_modules.o    \
          tglf_pkg.o        \
          tglf_allocate.o   \
          tglf_deallocate.o \
          tglf_startup.o    \
          tglf_hermite.o    \
          tglf_inout.o      \
          tglf_setup_geometry.o \
          tglf_LS.o         \
          tglf_eigensolver.o \
          tglf_geometry.o   \
          tglf_matrix.o     \
          tglf_max.o        \
          tglf_interface.o  \
          tglf_error.o      \
          tglf_isinf.o      \
          tglf_shutdown.o   \
          tglf_read_input.o \
          tglf_multiscale_spectrum.o \
          tglf_kygrid.o     \
	  tglf_run.o        \
          tglf_TM.o         \
          tglf_init_mpi.o   \
          tglf_run_mpi.o    \
          tglf_TM_mpi.o     \
          tglf_harvest.o    \
          tglf_nn.o         \
          tglf_run_nn.o     \
          tglf_run_nn_mpi.o \

.SUFFIXES : .o .f90 .F90

all: $(LLIB).a $(EXEC1) $(EXEC2) $(EXEC3)

$(EXEC1): $(LLIB).a $(EXEC1).o $(EXTRA_LIBS)
	$(FC) $(FFLAGS) -o $(EXEC1) $(EXEC1).o $(LLIB).a $(EXTRA_LIBS) $(LMATH)

$(EXEC2): $(LLIB).a $(EXEC2).o $(EXTRA_LIBS)
	$(FC) $(FFLAGS) -o $(EXEC2) $(EXEC2).o $(LLIB).a $(EXTRA_LIBS) $(LMATH)

$(EXEC3): $(LLIB).a $(EXEC3).o $(EXTRA_LIBS)
	$(FC) $(FFLAGS) -o $(EXEC3) $(EXEC3).o $(LLIB).a $(EXTRA_LIBS) $(LMATH)

$(LLIB).a: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

tglf_run_nn_mpi.o : tglf_run_nn.F90
	MPI_TGLF=1; $(FC) -o tglf_run_nn_mpi.o $(FMATH) $(FFLAGS) -c $<

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

.F90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o *.lst $(LLIB).a $(EXEC1) $(EXEC2) $(EXEC3)
	rm -f ${GACODE_ROOT}/modules/tglf*.mod 
