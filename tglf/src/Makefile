include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

ifeq ($(MPI),true)
   EXEC1=xtglf_mpi
   EXEC2=tglf_mpi
   SRC1=tglf_driver_mpi
   SRC2=tglf_mpi
else
   EXEC1=xtglf
   EXEC2=tglf
   SRC1=tglf_driver
   SRC2=tglf
endif

LLIB=tglf_lib

OBJECTS = tglf_modules.o   \
          tglf_pkg.o       \
          tglf_LS.o        \
          tglf_geometry.o  \
          tglf_matrix.o    \
          tglf_max.o       \
          tglf_interface.o \
          tglf_init.o      \
          tglf_run.o       \
          tglf_run_mpi.o   \
          tglf_read_input.o 

.SUFFIXES : .o .f90

all: $(LLIB).a $(EXEC1) $(EXEC2)

$(EXEC1): $(LLIB).a $(SRC1).o
	$(FC) $(FFLAGS) -o $(EXEC1) $(SRC1).o $(LLIB).a $(EXTRA_LIBS)

$(EXEC2): $(LLIB).a $(SRC2).o
	$(FC) $(FFLAGS) -o $(EXEC2) $(SRC2).o $(LLIB).a $(EXTRA_LIBS)

$(LLIB).a: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(LLIB).a xtglf xtglf_mpi tglf tglf_mpi
	rm -f ${GACODE_ROOT}/modules/tglf*.mod 
