include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

OBJECTS =  hdf5_fcapi.o 

.SUFFIXES : .o .f90  .F90 

#Always create dummy and then overwite if needed
all: vshdf5_dummy $(VSHDF5_MAKECMD)

vshdf5: $(OBJECTS)
	$(ARCH) vshdf5_lib.a $(OBJECTS)

vshdf5_dummy: dummy.o
	$(ARCH) vshdf5_lib.a dummy.o

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) $(HDF5_INC) -c $<

.F90.o :
	$(FC) $(FMATH) $(FFLAGS) $(HDF5_INC) -c $<

#SEK: Need to figure this out.
clean:
	rm -f *.o $(LLIB).a ${GACODE_ROOT}/modules/expro_interface.mod
