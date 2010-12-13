include ../shared/install/make.inc.${GACODE_PLATFORM}
include install/make.ext.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG} ${GATOFLAGS}
else
   FFLAGS=${FOPT} ${GATOFLAGS}
endif


EXEC = gato

OBJECTS = $(EXEC).o gato_plotg_library.stub.o

.SUFFIXES : .o .f90

all: $(OBJECTS)
	$(FC) $(FMATH) $(FFLAGS) -o $(EXEC) $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

.f.o :
	$(F77) $(FFLAGS) -c $<

clean:
	rm -f *.o *~ $(EXEC)

include install/make.ext.${GACODE_PLATFORM}
