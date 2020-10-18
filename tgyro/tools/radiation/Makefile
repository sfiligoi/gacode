include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC=curves

OBJECTS = $(EXEC).o

.SUFFIXES : .o .f90

all: lib $(EXEC).o
	$(FC) $(FFLAGS) -o $(EXEC) $(EXEC).o ${GACODE_ROOT}/tgyro/src/tgyro_rad.o

lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -rf *.o $(EXEC)
