include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

export EXTRA_LIBS = \
	${NETCDF} \
	${GACODE_ROOT}/shared/fluxfit/fluxfit_lib.a \
   	${GACODE_ROOT}/shared/math/math_lib.a \
	${LMATH} 

EXEC = profiles_3d

OBJECTS = profiles_3d.o

.SUFFIXES : .o .f90

all: $(EXEC)

$(EXEC): $(OBJECTS)
	$(FC) $(FMATH) $(FFLAGS) -o $(EXEC) $(OBJECTS) $(EXTRA_LIBS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -I $(NETCDF_INC) -c $<

clean:
	rm -f *.o $(EXEC)



