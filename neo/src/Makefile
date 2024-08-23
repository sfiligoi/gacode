include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG} ${FOMP}
   CFLAGS = -g ${FOMP}
else
   FFLAGS=${FOPT} ${FOMP}
   CFLAGS = -O3 ${FOMP}
endif

EXEC = neo

ifdef FANN_ROOT
   # neural net libraries
   NN=neo_neural.fann
else
   #stub file
   NN=neo_neural
endif

OBJECTS = neo_globals.o \
          neo_energy_grid.o \
          neo_interface.o \
          neo_allocate_profile.o \
          neo_umfpack.o \
          neo_sparse_solve.o \
          neo_equilibrium.o \
          neo_g_velocitygrids.o \
          neo_rotation.o \
          neo_nclass_dr.o \
          neo_theory.o \
          neo_transport.o \
          neo_3d_driver.o \
          neo_check.o \
          neo_compute_fcoll.o \
	  $(NN).o \
          neo_do.o \
          neo_error.o \
          neo_init.o \
          neo_init_serial.o \
          neo_make_profiles.o \
          neo_read_input.o \
          neo_run.o \
          neo_spitzer.o \
          matconv.o

# JC: Commenting out 2024.04.25
#NEO_FLAGS = -cpp
#ifdef CUDA_DIR
#   OBJECTS   += neo_cusolve.o cusolve.o
#   NEO_FLAGS += -DNEO_HAVE_CUDA $(CUDA_INC)
#   NEO_LIBS  += $(CUDA_LIB)
#endif
#ifdef PETSC_DIR
#   OBJECTS   += neo_petsc.o
#   NEO_FLAGS += -DNEO_HAVE_PETSC $(PETSC_INC)
#   NEO_LIBS  += $(PETSC_LIB)
#endif
#ifdef SUPERLU_DIR
#   OBJECTS   += neo_superlu.o st_to_cc.o c_fortran_dgssv.o
#   NEO_FLAGS += -DNEO_HAVE_SUPERLU $(SUPERLU_INC)
#   NEO_LIBS  += $(SUPERLU_LIB)
#endif

.SUFFIXES : .o .f90 .F90

all: neo_lib.a $(EXEC).o $(EXTRA_LIBS)
	$(FC) $(FFLAGS) $(NEO_FLAGS) -o $(EXEC) $(EXEC).o neo_lib.a $(EXTRA_LIBS) $(LMATH) $(NN_LIB) ${NEO_LIBS}

neo_lib.a: $(OBJECTS)
	$(ARCH) neo_lib.a $(OBJECTS)

.f90.o .F90.o :
	$(FC) $(FMATH) $(FFLAGS) $(NN_LIB) $(NEO_FLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) $(CDEFS) $(NEO_FLAGS) -c $< $(VERBOSE)

clean:
	rm -rf *.o neo_lib.a $(EXEC) *~
	cd ${GACODE_ROOT}/modules ; rm -f neo*.mod matconv.mod

-include ../install/make.ext.${GACODE_PLATFORM}

# all module induced dependencies in directory /home/ksh/gacode/neo/src:

# module matconv in matconv.f90
neo_cusolve.o: matconv.o

# module neo_3d_driver in neo_3d_driver.f90
neo_do.o: neo_3d_driver.o

# module neo_allocate_profile in neo_allocate_profile.f90
neo_3d_driver.o neo_do.o neo_make_profiles.o: neo_allocate_profile.o

# module neo_energy_grid in neo_energy_grid.f90
neo_3d_driver.o neo_do.o neo_g_velocitygrids.o neo_spitzer.o neo_transport.o: neo_energy_grid.o

# module neo_equilibrium in neo_equilibrium.f90
neo_do.o neo_nclass_dr.o neo_neural.o neo_rotation.o neo_theory.o neo_transport.o: neo_equilibrium.o

# module neo_g_velocitygrids in neo_g_velocitygrids.f90
neo_3d_driver.o neo_do.o: neo_g_velocitygrids.o

# module neo_globals in neo_globals.f90
neo.o neo_3d_driver.o neo_allocate_profile.o neo_check.o neo_do.o neo_energy_grid.o neo_equilibrium.o neo_error.o neo_g_velocitygrids.o neo_init.o neo_init_serial.o neo_interface.o neo_make_profiles.o neo_nclass_dr.o neo_neural.o neo_petsc.o neo_rbf.o neo_read_input.o neo_rotation.o neo_run.o neo_sparse_solve.o neo_spitzer.o neo_superlu.o neo_theory.o neo_transport.o neo_umfpack.o: neo_globals.o

# module neo_interface in neo_interface.f90
neo_init.o neo_init_serial.o neo_run.o: neo_interface.o

# module neo_nclass_dr in neo_nclass_dr.f90
neo_do.o neo_theory.o: neo_nclass_dr.o

# module neo_neural in $(NN).f90
neo_do.o: $(NN).o

# module neo_rotation in neo_rotation.f90
neo_do.o neo_transport.o: neo_rotation.o

# module neo_sparse_solve in neo_sparse_solve.F90
neo_do.o: neo_sparse_solve.o

# module neo_theory in neo_theory.f90
neo_do.o neo_transport.o: neo_theory.o

# module neo_transport in neo_transport.f90
neo_do.o: neo_transport.o

# module neo_umfpack in neo_umfpack.F90
neo_sparse_solve.o neo_spitzer.o: neo_umfpack.o
