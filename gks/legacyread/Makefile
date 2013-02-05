include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

LLIB=legacyread_lib

OBJECTS = legacyread_interface.o \
          legacyread_avg.o      \
          legacyread_run.o       \
          legacyread_ufiles.o \
          legacyread_iterdb.o   \
	  legacyread_gridsetup.o  \
          legacyread_LIBPACK.o 

.SUFFIXES : .o .f .f90

lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.lst $(LLIB).a 
	rm -f ${GACODE_ROOT}/modules/legacyread*.mod
