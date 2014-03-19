include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif
      
EXEC = vgen

OBJECTS = vgen_globals.o \
	  vgen_init.o \
	  vgen_compute_neo.o \
	  vgen_getgeo.o \
          vgen_reduce.o \
	  $(EXEC).o 

.SUFFIXES : .o .f90

all: $(OBJECTS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS) ../../src/neo_lib.a $(EXTRA_LIBS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(EXEC) ${GACODE_ROOT}/modules/vgen_globals.mod
