# Makefile for gyro/vshdf5
include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

OBJECTS = \
	gyro_vshdf5_mod.o \
	gyro_write_master_hdf5.o 

DUM_OBJECTS = \
	gyro_dummy_master_hdf5.o

.SUFFIXES : .o .f90 

# Always make the dummy.  If we really want the library it will be
# overwritten
all: gyro_vshdf5_dummy $(VSHDF5_MAKECMD) 


VSHDF5_INC=-I../../shared/vshdf5

vshdf5: $(OBJECTS)
	rm gyrovshdf5_lib.a
	$(ARCH) gyrovshdf5_lib.a $(OBJECTS)

gyro_vshdf5_dummy: $(DUM_OBJECTS)
	$(ARCH) gyrovshdf5_lib.a $(DUM_OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) $(VSHDF5_INC) $(HDF5_INC) $(VSHDF5_INC) -c $<

clean:  
	rm -f *.o *.a

include ${GACODE_ROOT}/gyro/install/make.ext.${GACODE_PLATFORM}
