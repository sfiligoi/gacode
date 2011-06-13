include ../../shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC = neo

OBJECTS = neo_globals.o \
          neo_error.o \
	  neo_profile_exp.o \
          neo_check.o \
          neo_sparse_solve.o \
          neo_compute_fcoll.o \
          neo_energy_grid.o \
	  neo_allocate_profile.o \
	  neo_make_profiles.o \
          neo_read_input.o \
          neo_equilibrium.o \
	  neo_rotation.o \
	  neo_nclass_dr.o \
	  neo_theory.o \
	  neo_transport.o \
          neo_experimental_profiles.o \
          neo_map_experimental_profiles.o \
	  neo_g_velocitygrids.o \
          neo_interface.o \
          neo_init.o \
          neo_run.o \
          neo_spitzer.o \
          neo_do.o

.SUFFIXES : .o .f90

all: neo_lib.a $(EXEC).o
	$(FC) $(FFLAGS) -o $(EXEC) $(EXEC).o neo_lib.a $(EXTRA_LIBS)

neo_lib.a: $(OBJECTS)
	$(ARCH) neo_lib.a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:  
	rm -rf *.o neo_lib.a $(EXEC) *~
