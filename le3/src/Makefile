include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC = le3

OBJECTS = le3_globals.o \
          le3_func.o \
	  le3_alloc.o \
          le3_basis.o \
          le3_compute_theory.o \
          le3_map.o \
          le3_rz.o \
          le3_geometry.o \
          le3_geometry_rho.o \
	  le3_geometry_matrix.o \
          $(EXEC).o

.SUFFIXES : .o .f90

all: $(OBJECTS) $(EXTRA_LIBS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS) $(EXTRA_LIBS) $(LMATH)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(EXEC) 

# all module induced dependencies in directory /home/ksh/gacode/le3/src:

# module le3_globals in le3_globals.f90
le3.o le3_alloc.o le3_basis.o le3_compute_theory.o le3_func.o le3_geometry.o le3_geometry_matrix.o le3_geometry_rho.o le3_rz.o: le3_globals.o
