#-----------------------------------------------------
# Makefile to compile the Dep_driver test system.
# Written to work on NERSC Hopper machine
#-----------------------------------------------------

# Define compilers and flags.
include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}
#FC      = ftn
#FOPT    = -fast
#FDEBUG  = -g -Ktrap=fp -Mbounds -Minform,inform
export EXTRA_LIBS = ${GACODE_ROOT}/tglf/src/tglf_lib.a

#ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
#else
#   FFLAGS=${FOPT}
#endif

EXEC=Dep_driver

OBJECTS = from_tglf_to_Dep.o \
          from_nubeam_to_Dep.o \
          from_Dep_to_nubeam.o \
          Dep_global.o \
          Dep_grid_wts.o \
          Dep_kernel.o \
          Dep_read_input.o \
          i0.o \
          j0y0.o \
          Dep_allocate.o \
          Dep_mainsub.o \
          Dep_driver.o

.SUFFIXES : .o .f90 .f

all: $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) $(EXTRA_LIBS) -o $(EXEC)

.f90.o :
	$(FC) $(FFLAGS) -r8 -c $<

.f.o :
	$(FC) $(FFLAGS) -r4 -c $<

clean:
	rm -f *.o *.mod $(EXEC)
