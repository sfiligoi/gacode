# Makefile for gyro/vshdf5
include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

LLIB=gyrovshdf5_lib

ifeq ($(VSHDF5_MAKECMD),vshdf5)
   OBJECTS = \
	gyro_vshdf5_mod.o \
	gyro_write_initdata_hdf5.o \
	gyro_write_timedata_hdf5.o \
	gyro_write_timedata_wedge_hdf5.o
else
   OBJECTS = \
	gyro_write_dummy_hdf5.o
endif

.SUFFIXES : .o .f90 

VSHDF5_INC=-I../../shared/vshdf5

lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) $(VSHDF5_INC) $(HDF5_INC) $(VSHDF5_INC) -c $<

clean:  
	rm -f *.o *.a ${GACODE_ROOT}/modules/*vshdf5*.mod

include ${GACODE_ROOT}/gyro/install/make.ext.${GACODE_PLATFORM}
