include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC1=tglf
EXEC2=tglf_mpi

LLIB=tglf_lib

ifdef FANN_ROOT
   # neural-net file
   NN=tglf_nn_TM.fann
else
   # stub file
   NN=tglf_nn_TM
endif

OBJECTS = tglf_isnan.o      \
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
          $(NN).o

.SUFFIXES : .o .f90 .F90

all: $(LLIB).a $(EXEC1) $(EXEC2) Makefile

$(EXEC1): $(LLIB).a $(EXEC1).o
	$(FC) $(FFLAGS) -o $(EXEC1) $(EXEC1).o $(LLIB).a $(LMATH) $(NN_LIB)
$(EXEC2): $(LLIB).a $(EXEC2).o
	$(FC) $(FFLAGS) -o $(EXEC2) $(EXEC2).o $(LLIB).a $(LMATH) $(NN_LIB)

$(LLIB).a: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

tglf_run_mpi.o: tglf_run.F90
	$(FC) -o tglf_run_mpi.o $(FMATH) $(FFLAGS) -DMPI_TGLF=1 -c $<

.f90.o:
	$(FC) $(FMATH) $(FFLAGS) $(NN_LIB) -c $<
.F90.o:
	$(FC) $(FMATH) $(FFLAGS) $(NN_LIB) -c $<

clean:
	rm -f *.o *.lst $(LLIB).a $(EXEC1) $(EXEC2)
	rm -f ${GACODE_ROOT}/modules/tglf*.mod

-include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

# all module induced dependencies in directory tglf/src

# module tglf_tg in tglf_TM_driver.f90
tglf_TM_driver.o tglf_driver.o tglf_driver_mpi.o: tglf_TM_driver.o

# module tglf_mpi in tglf_init_mpi.f90
tglf_TM_mpi.o tglf_driver_mpi.o tglf_init_mpi.o tglf_mpi.o tglf_nn_TM.o tglf_run.o tglf_run_mpi.o: tglf_init_mpi.o

# module tglf_interface in tglf_interface.f90
tglf.o tglf_init_mpi.o tglf_mpi.o tglf_nn_TM.o tglf_read_input.o tglf_run.o tglf_run_mpi.o: tglf_interface.o

# module tglf_max_dimensions in tglf_modules.f90
tglf_interface.o tglf_modules.o: tglf_modules.o

# module tglf_dimensions in tglf_modules.f90
tglf_LS.o tglf_TM.o tglf_TM_mpi.o tglf_allocate.o tglf_eigensolver.o tglf_geometry.o tglf_hermite.o tglf_inout.o tglf_kygrid.o tglf_matrix.o tglf_max.o tglf_modules.o tglf_multiscale_spectrum.o tglf_nn_TM.o tglf_pkg.o tglf_setup_geometry.o tglf_startup.o: tglf_modules.o

# module tglf_global in tglf_modules.f90
tglf.o tglf_LS.o tglf_TM.o tglf_TM_mpi.o tglf_eigensolver.o tglf_geometry.o tglf_hermite.o tglf_inout.o tglf_interface.o tglf_kygrid.o tglf_matrix.o tglf_max.o tglf_mpi.o tglf_multiscale_spectrum.o tglf_nn_TM.o tglf_pkg.o tglf_run.o tglf_run_mpi.o tglf_setup_geometry.o tglf_startup.o: tglf_modules.o

# module tglf_closure in tglf_modules.f90
tglf_eigensolver.o: tglf_modules.o

# module tglf_hermite in tglf_modules.f90
tglf_geometry.o tglf_hermite.o tglf_matrix.o: tglf_modules.o

# module tglf_species in tglf_modules.f90
tglf_LS.o tglf_TM.o tglf_TM_mpi.o tglf_allocate.o tglf_deallocate.o tglf_eigensolver.o tglf_geometry.o tglf_inout.o tglf_kygrid.o tglf_matrix.o tglf_max.o tglf_multiscale_spectrum.o tglf_setup_geometry.o tglf_startup.o: tglf_modules.o

# module tglf_kyspectrum in tglf_modules.f90
tglf_TM.o tglf_TM_mpi.o tglf_inout.o tglf_kygrid.o tglf_multiscale_spectrum.o: tglf_modules.o

# module tglf_eigen in tglf_modules.f90
tglf_LS.o tglf_eigensolver.o: tglf_modules.o

# module tglf_weight in tglf_modules.f90
tglf_LS.o: tglf_modules.o

# module tglf_xgrid in tglf_modules.f90
tglf_TM.o tglf_TM_mpi.o tglf_geometry.o tglf_matrix.o tglf_multiscale_spectrum.o: tglf_modules.o

# module tglf_sgrid in tglf_modules.f90
tglf_LS.o tglf_eigensolver.o tglf_geometry.o tglf_inout.o tglf_matrix.o: tglf_modules.o

# module tglf_coeff in tglf_modules.f90
tglf_LS.o tglf_allocate.o tglf_deallocate.o tglf_eigensolver.o tglf_inout.o tglf_matrix.o: tglf_modules.o

# module tglf_nuei_coeff in tglf_modules.f90
tglf_eigensolver.o: tglf_modules.o

# module tglf_tg in tglf_modules.f90
tglf_TM_driver.o tglf_driver.o tglf_driver_mpi.o: tglf_modules.o

# module tglf_pkg in tglf_pkg.f90
tglf.o tglf_TM_driver.o tglf_driver.o tglf_driver_mpi.o tglf_inout.o tglf_max.o tglf_mpi.o tglf_nn_TM.o tglf_run.o tglf_run_mpi.o: tglf_pkg.o
