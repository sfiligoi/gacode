include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif
      
EXEC = squeeze

OBJECTS = $(EXEC).o 

.SUFFIXES : .o .f90

all: $(OBJECTS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS) ../../src/neo_lib.a $(EXTRA_LIBS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(EXEC)
