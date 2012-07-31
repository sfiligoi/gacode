#-----------------------------------------------------
# Makefile to compile the Dep_driver test system.
#-----------------------------------------------------

# Define compilers and flags.
include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}
export EXTRA_LIBS = ${GACODE_ROOT}/tglf/src/tglf_lib.a

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC=Dep_driver

LLIB=dep_lib

OBJECTS = Dep_from_tglf.o \
          Dep_from_nubeam.o \
          Dep_to_nubeam.o \
          Dep_global.o \
          Dep_grid_wts.o \
          Dep_kernel.o \
          Dep_read_input.o \
          Dep_i0.o \
          Dep_j0y0.o \
          Dep_allocate.o \
          Dep_mainsub.o 

.SUFFIXES : .o .f90 .f

all: $(LLIB).a $(EXEC)

$(EXEC): $(LLIB).a $(EXEC).o
	 $(FC) $(FFLAGS) -o $(EXEC) $(EXEC).o $(LLIB).a $(EXTRA_LIBS)

$(LLIB).a: $(OBJECTS)
	 $(ARCH) $(LLIB).a $(OBJECTS) 

.f90.o :
	$(FC) $(FFLAGS) -r8 -c $<

.f.o :
	$(FC) $(FFLAGS) -r4 -c $<

clean:
	rm -f *.o  $(EXEC) $(LLIB).a
	rm -f ${GACODE_ROOT}/modules/Dep*.mod
