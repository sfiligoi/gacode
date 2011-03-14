include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

OBJECTS =  hdf5_fcapi.o 

.SUFFIXES : .o .f90  .F90 

all: $(VSHDF5_MAKECMD)

vshdf5: $(OBJECTS)
	$(ARCH) libvshdf5.a $(OBJECTS)

.f.o :
	$(FC) $(FMATH) $(FFLAGS) $(HDF5_INC) -c $<

.F.o :
	$(FC) $(FMATH) $(FFLAGS) $(HDF5_INC) -c $<

#SEK: Need to figure this out.
clean:
	rm -f *.o $(LLIB).a ${GACODE_ROOT}/modules/expro_interface.mod
