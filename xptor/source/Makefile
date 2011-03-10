include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

LLIB=../driver/xptor_lib

OBJECTS = neucg.o \
          neuden.o \
          neut.o \
	  palpha_dv.o \
	  prad_dv.o \
	  pohmic_dv.o \
          sigvth.o  

.SUFFIXES : .o .f 

lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS) 

.f.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o
