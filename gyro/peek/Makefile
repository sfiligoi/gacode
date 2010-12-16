include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC = peek

OBJ = peek_globals.o \
      read_profile_vugyro.o \
      read_time.o \
      $(EXEC).o

.SUFFIXES : .o .f90 

all: $(OBJ)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJ) ../src/gyro_lib.a

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:  
	rm -f *.o $(EXEC)
	rm -f ${GACODE_ROOT}/modules/peek_globals.mod
