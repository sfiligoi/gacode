include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC=reader
LLIB=dataread_lib

OBJECTS = data_interface.o \
          data_run.o       \
          data_readufiles.o \
          data_readiterdb.o   \
          data_avg.o      \
	  data_gridsetup.o  \
          data_LIBPACK.o 

.SUFFIXES : .o .f .f90

all: $(LLIB).a $(EXEC) prgen_read_ufile.o

$(LLIB).a: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

$(EXEC): prgen_read_ufile.o
	$(FC) $(FMATH) -o $(EXEC) prgen_read_ufile.o $(LLIB).a

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o $(LLIB).a 
	rm -f ${GACODE_ROOT}/modules/data_interface.mod
