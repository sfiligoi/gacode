include ${GYRO_DIR}/install/make.inc.${GYRO_PLAT}

ifeq ($(MODE),DEBUG)
   FFLAGS=$(FMATH) $(FDEBUG)
else
   FFLAGS=$(FMATH) $(FOPT)
endif

all: spec_xy spec_n spec_i

OBJECTS = input_data.o \
          read_input.o \
          per_deriv.o

.SUFFIXES : .o .f90

.f90.o :
	$(FC) $(FFLAGS) -c $<

spec_xy: $(OBJECTS)
	$(FC) $(FFLAGS) -o spec_xy spec_xy.f90 $(OBJECTS)

spec_n: $(OBJECTS)
	$(FC) $(FFLAGS) -o spec_n spec_n.f90 $(OBJECTS)

spec_i: $(OBJECTS)
	$(FC) $(FFLAGS) -o spec_i spec_i.f90 $(OBJECTS)

clean:
	rm -f *.o spec_xy spec_n spec_i *~

# all module induced dependencies in directory gyro/tools/spectra/src

# module input_data in input_data.f90
read_input.o spec_i.o spec_n.o spec_xy.o: input_data.o
