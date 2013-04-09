include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG} 
else
   FFLAGS=${FOPT} 
endif

ifeq ($(HDF5),true)
  FFLAGS += ${FHDF5} ${HDF5_FLAGS}
  OBJECTS = hdf5_fcapi.o 
else
  OBJECTS = hdf5_fcdummy.o 
endif

LLIB = vshdf5_lib



.SUFFIXES : .o .f90  .F90 

##Always create dummy and then overwite if needed
#all: vshdf5_dummy $(VSHDF5_MAKECMD)
#
#vshdf5: $(OBJECTS)
#	$(ARCH) vshdf5_lib.a $(OBJECTS)
#
#vshdf5_dummy: dummy.o
#	$(ARCH) vshdf5_lib.a dummy.o
#
lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) $(HDF5_INC) -c $<

.F90.o :
	$(FC) $(FMATH) $(FFLAGS) $(HDF5_INC) -c $<

clean:
	rm -f *.o *.a ${GACODE_ROOT}/modules/*hdf5*.mod
