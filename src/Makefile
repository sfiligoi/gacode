include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC = cgyro

OBJECTS = cgyro_globals.o \
          cgyro_io.o \
          cgyro_mpi_grid.o \
          cgyro_check.o \
	  cgyro_make_profiles.o \
          cgyro_read_input.o \
          cgyro_equilibrium.o \
          cgyro_interface.o \
          cgyro_init.o \
          cgyro_run.o \
          cgyro_gyro.o \
          cgyro_field.o \
          cgyro_freq.o \
          cgyro_collision.o \
	  cgyro_gk.o \
          cgyro_do.o \
          cgyro_write_timedata.o 

.SUFFIXES : .o .f90

all: cgyro_lib.a $(EXEC).o
	$(FC) $(FFLAGS) -o $(EXEC) $(EXEC).o cgyro_lib.a $(EXTRA_LIBS)

cgyro_lib.a: $(OBJECTS)
	$(ARCH) cgyro_lib.a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:  
	rm -rf *.o cgyro_lib.a $(EXEC) *~
