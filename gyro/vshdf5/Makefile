# Makefile for gyro/vshdf5
include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

OBJECTS = \
	gyro_hdf5_mod.f90 \
	gyro_write_master_hdf5.f90

DUM_OBJECTS = \
	write_hdf5_dummy.f90

.SUFFIXES : .o .f90 

all: $(VSHDF5_MAKECMD) gyro_vshdf5_dummy

VSHDF5_INC=-I../../shared/vshdf5

vshdf5: $(OBJECTS)
	$(ARCH) libgyro_vshdf5.a $(OBJECTS)

gyro_vshdf5_dummy: $(OBJECTS)
	$(ARCH) libgyro_vshdf5_dummy.a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) $(VSHDF5_INC) -c $<

clean:  
	rm -f *.o $(LLIB).a
