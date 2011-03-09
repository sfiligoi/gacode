include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

LLIB=../driver/xptor_lib

OBJECTS = advancefields.o \
          norm_dv.o \
          perturb_dv.o \
          ptor_dv.o \
          ptorinit_dv.o \
          trcoef_dv.o \
          vtrans_dv.o \
          tglf_startup.o

.SUFFIXES : .o .f 


lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS) 

.f.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o 
