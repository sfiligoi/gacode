include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

LLIB=qfm_lib

OBJECTS = qfm_sub.o \
          qfm_interface.o \
          qfm_run.o \
          difftd.o \
          disp9t.o \
          r8tomsqz.o

.SUFFIXES : .o .f .f90

all: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

.f.o :
	$(F77) $(FFLAGS) -c $<

clean:
	rm -f *.o $(LLIB).a ${GACODE_ROOT}/modules/qfm*.mod
