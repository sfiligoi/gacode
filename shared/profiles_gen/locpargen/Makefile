include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC = locpargen

OBJECTS = locpargen_globals.o \
          $(EXEC).o

.SUFFIXES : .o .f90

exec: $(OBJECTS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(EXEC)
	rm -f ${GACODE_ROOT}/modules/locpargen*.mod
