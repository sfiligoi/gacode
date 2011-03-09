include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC1=xtglf
EXEC2=tglf

LLIB=tglf_lib

OBJECTS = tglf_modules.o \
          tglf_LS.o       \
          tglf_geometry.o \
          tglf_matrix.o   \
          tglf_max.o      \
          tglf_pkg.o \
          tglf_interface.o \
          tglf_init.o \
          tglf_run.o \
          tglf_read_input.o

.SUFFIXES : .o .f90

all: $(LLIB).a $(EXEC1) $(EXEC2)

$(EXEC1): $(LLIB).a tglf_driver.o
	$(FC) $(FFLAGS) -o $(EXEC1) tglf_driver.o $(LLIB).a $(EXTRA_LIBS)

$(EXEC2): $(LLIB).a tglf.o
	$(FC) $(FFLAGS) -o $(EXEC2) tglf.o $(LLIB).a $(EXTRA_LIBS)

$(LLIB).a: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(LLIB).a $(EXEC1) $(EXEC2)
	rm -f ${GACODE_ROOT}/modules/tglf*.mod 
