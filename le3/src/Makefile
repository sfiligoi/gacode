include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC = le3

OBJECTS = le3_globals.o \
          le3_func.o \
	  le3_alloc.o \
          le3_map.o \
          le3_rz.o \
          le3_write.o \
          $(EXEC).o

.SUFFIXES : .o .f90

all: $(OBJECTS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS) $(EXTRA_LIBS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(EXEC) 
