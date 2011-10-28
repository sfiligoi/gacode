include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else 
   FFLAGS=${FOPT}
endif      

LLIB=GKEIGEN_lib

OBJECTS = GKEIGEN_globals.o \
          GKEIGEN_do.o \
          GKEIGEN_matrix_read.o \
          GKEIGEN_matrix_write.o \
          GKEIGEN_catch_error.o \
          GKEIGEN_split_comm.o

ifeq ($(OMP),true)
   FC += ${FOMP}
else
   FC += ${FNOMP}
   OBJECTS += ../src/gyro_omp_lib.o
endif
                  
.SUFFIXES : .o .f90 .F

lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -rf *.o $(LLIB).a
	rm -rf ${GACODE_ROOT}/modules/gkeigen_globals.mod	                

include ${GACODE_ROOT}/gyro/install/make.ext.${GACODE_PLATFORM}
