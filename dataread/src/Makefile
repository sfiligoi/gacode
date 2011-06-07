include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif


LLIB=dataread_lib

OBJECTS = data_interface.o \
          data_run.o       \
          data_readufiles.o \
          data_readiterdb.o   \
          data_avg.o      \
	  data_gridsetup.o  \
          data_LIBPACK.o

.SUFFIXES : .o .f .f90

all: $(LLIB).a 

$(LLIB).a: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<
.f.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(LLIB).a 
	rm -f ${GACODE_ROOT}/modules/data_interface.mod
