include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC=xtglf
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

all: $(LLIB).a tglf_driver.o
	$(FC) $(FFLAGS) -o $(EXEC) tglf_driver.o $(LLIB).a $(EXTRA_LIBS)

tglf_lib.a: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(LLIB).a $(EXEC)
	rm -f ${GACODE_ROOT}/modules/tglf*.mod 
