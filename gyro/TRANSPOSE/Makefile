include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

FC += ${FOMP}

LLIB=TRANSP_lib

OBJECTS = fTRANSP_GLOBALS.o \
          fTRANSP_INIT.o \
          fTRANSP_DO.o \
          fTRANSP_CLEANUP.o \
          rTRANSP_GLOBALS.o \
          rTRANSP_INIT.o \
          rTRANSP_DO.o \
          rTRANSP_CLEANUP.o 

.SUFFIXES : .o .f90

lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS) 

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:  
	rm -rf *.o $(LLIB).a
	rm -rf ${GACODE_ROOT}/modules/ftransp*.mod
	rm -rf ${GACODE_ROOT}/modules/rtransp*.mod

# all module induced dependencies in directory gyro/TRANSPOSE

# module fTRANSP_GLOBALS in fTRANSP_GLOBALS.f90
fTRANSP_CLEANUP.o fTRANSP_DO.o fTRANSP_INIT.o: fTRANSP_GLOBALS.o

# module rTRANSP_GLOBALS in rTRANSP_GLOBALS.f90
rTRANSP_CLEANUP.o rTRANSP_DO.o rTRANSP_INIT.o: rTRANSP_GLOBALS.o
