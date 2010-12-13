include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC=tglf

OBJECTS = tglf_modules.o \
          tglf_pkg.o \
          tglf_interface.o \
          tglf_init.o \
          tglf_run.o \
          tglf_read_input.o

.SUFFIXES : .o .f90

      all: tglf_lib.a $(EXEC).o
	$(FC) $(FFLAGS) -o $(EXEC) $(EXEC).o tglf_lib.a $(EXTRA_LIBS)

tglf_lib.a: $(OBJECTS)
	$(ARCH) tglf_lib.a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -rf *.o tglf_lib.a $(EXEC) *~
