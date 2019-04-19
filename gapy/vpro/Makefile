include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

FC += ${FOMP}

OBJECTS = gapy_vpro.o \
          gapy_vpro_util.o \
          gapy_vpro_comm.o \
          expro_locsim_interface.o \
          expro_locsim_profiles.o \
          expro_locsim_alloc.o 


.SUFFIXES : .o .f90

all: vpro_lib.a

vpro_lib.a: $(OBJECTS)
	$(ARCH) vpro_lib.a $(OBJECTS)

.f90.o :
	$(FC) $(FFLAGS) $(FMATH) -c $<

clean:
	rm -f *.o *.a *.mod
