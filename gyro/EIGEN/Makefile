include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else 
   FFLAGS=${FOPT}
endif      

LLIB=EIGEN_lib

OBJECTS = EIGEN_globals.o \
          EIGEN_do.o \
          EIGEN_matrix_read.o \
          EIGEN_matrix_write.o

.SUFFIXES : .o .f90 .F

lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -rf *.o $(LLIB).a
	rm -rf ${GACODE_ROOT}/modules/eigen_globals.mod	                

include ${GACODE_ROOT}/gyro/install/make.ext.${GACODE_PLATFORM}
