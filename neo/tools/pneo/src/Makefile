include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC1 = pneo
EXEC2 = pneo_rbf

OBJECTS = pneo_globals.o \
          $(EXEC1).o

.SUFFIXES : .o .f90

all: exec1 exec2

exec1: $(OBJECTS)
	$(FC) $(FFLAGS) -o $(EXEC1) $(OBJECTS) ${GACODE_ROOT}/neo/src/neo_lib.a $(EXTRA_LIBS) $(LMATH) $(NN_LIB) 


exec2: $(OBJECTS) $(EXEC2).o 
	$(FC) $(FFLAGS) -o $(EXEC2) $(EXEC2).o $(LMATH)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(EXEC1) $(EXEC2) ${GACODE_ROOT}/modules/pneo*.mod
