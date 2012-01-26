include ../../shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC = gkcoll

OBJECTS = gkcoll_globals.o \
          gkcoll_error.o \
	  gkcoll_profile_exp.o \
          gkcoll_check.o \
	  gkcoll_allocate_profile.o \
	  gkcoll_make_profiles.o \
          gkcoll_read_input.o \
          gkcoll_equilibrium.o \
          gkcoll_experimental_profiles.o \
          gkcoll_map_experimental_profiles.o \
          gkcoll_interface.o \
          gkcoll_init.o \
          gkcoll_init_serial.o \
          gkcoll_run.o \
	  gkcoll_legendre.o \
          gkcoll_gyro.o \
          gkcoll_poisson.o \
          gkcoll_freq.o \
          gkcoll_gk.o \
	  gkcoll_implicit.o \
	  gkcoll_allimplicit.o \
	  gkcoll_collision.o \
          gkcoll_do.o

.SUFFIXES : .o .f90

all: gkcoll_lib.a $(EXEC).o
	$(FC) $(FFLAGS) -o $(EXEC) $(EXEC).o gkcoll_lib.a $(EXTRA_LIBS)

gkcoll_lib.a: $(OBJECTS)
	$(ARCH) gkcoll_lib.a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:  
	rm -rf *.o gkcoll_lib.a $(EXEC) *~
