include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

LLIB=gyro_vshdf5
DUMLLIB=gyro_vshdf5_dummy

OBJECTS =  \
	hdf5_fcapi.o \
	hdf5_mod.o \
	write_distributed_complex_h5.o \
	write_distributed_real_h5.o \
	write_hdf5.o

DUM_OBJECTS =  \
	write_hdf5_dummy.o

.SUFFIXES : .o .f90  .F90 

all: $(VSHDF_MAKE)

gyro_vshdf5: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

gyro_vshdf5_dummy: $(DUM_OBJECTS)
	$(ARCH) $(DUMLLIB).a $(DUM_OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

.F90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

#SEK: Need to figure this out.
clean:
	rm -f *.o $(LLIB).a ${GACODE_ROOT}/modules/expro_interface.mod

