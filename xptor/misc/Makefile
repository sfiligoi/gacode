include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}
include ../install/make.ext.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG} ${XPTORFLAGS}
else
   FFLAGS=${FOPT} ${XPTORFLAGS}
endif

LLIB=../driver/xptor_lib

OBJECTS = r8tomsqz.o \
          inter_cspl.o \
          lspolysmooth.o \
          part.o \
          volint.o \
          volintk.o \
          pmisc.o   \
          pchmisc.o \
          curvss.o  \
          dpchfe.o \
          dpchim.o \
          dchfev.o \
          dgbtrs.o \
          dgbtrf.o \
          w_lin_interp_r8.o 

.SUFFIXES : .o .f .F


lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS) 

.f.o :
	$(F77) $(FFLAGS) -c $<

.F.o :
	$(F77) $(FFLAGS) -c $< f77_dcomplx.h

clean:
	rm -f *.o *.lst *.gch
