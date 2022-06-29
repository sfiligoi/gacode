include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC1=gftm
EXEC2=gftm_mpi

LLIB=gftm_lib

OBJECTS = gftm_isnan.o      \
          gftm_modules.o    \
          gftm_pkg.o        \
          gftm_allocate.o   \
          gftm_deallocate.o \
          gftm_startup.o    \
          gftm_hermite.o    \
          gftm_inout.o      \
          gftm_setup_geometry.o \
          gftm_LS.o         \
          gftm_eigensolver.o \
          gftm_geometry.o   \
          gftm_matrix.o     \
          gftm_max.o        \
          gftm_interface.o  \
          gftm_error.o      \
          gftm_isinf.o      \
          gftm_shutdown.o   \
          gftm_read_input.o \
          gftm_multiscale_spectrum.o \
          gftm_kygrid.o     \
          gftm_run.o        \
          gftm_TM.o         \
          gftm_init_mpi.o   \
          gftm_run_mpi.o    \
          gftm_TM_mpi.o     \
          gftm_collision1.o     \
          gftm_collision2.o


.SUFFIXES : .o .f90 .F90

all: $(LLIB).a $(EXEC1) $(EXEC2) Makefile

$(EXEC1): $(LLIB).a $(EXEC1).o
	$(FC) $(FFLAGS) -o $(EXEC1) $(EXEC1).o $(LLIB).a $(LMATH) 
$(EXEC2): $(LLIB).a $(EXEC2).o
	$(FC) $(FFLAGS) -o $(EXEC2) $(EXEC2).o $(LLIB).a $(LMATH) $(NN_LIB)

$(LLIB).a: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

gftm_run_mpi.o: gftm_run.F90
	$(FC) -o gftm_run_mpi.o $(FMATH) $(FFLAGS) -DMPI_GFTM=1 -c $<

.f90.o:
	$(FC) $(FMATH) $(FFLAGS) -c $<
.F90.o:
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o *.lst $(LLIB).a $(EXEC1) $(EXEC2)
	rm -f ${GACODE_ROOT}/modules/gftm*.mod

-include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}
