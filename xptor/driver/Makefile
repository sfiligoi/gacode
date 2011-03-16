include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

export EXTRA_LIBS =  \
       ${GACODE_ROOT}/tglf/src/tglf_lib.a  \
       ${GACODE_ROOT}/neo/src/neo_lib.a \
       ${GACODE_ROOT}/shared/EXPRO/EXPRO_lib.a \
       ${GACODE_ROOT}/shared/math/math_lib.a \
       ${GACODE_ROOT}/shared/GEO/GEO_lib.a \
       ${GACODE_ROOT}/shared/UMFPACK/UMFPACK_lib.a \
       ${GACODE_ROOT}/neo/nclass/nclass_lib.a \
       ${LMATH} \
       ${NETCDF}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC= xptor
LLIB=xptor_lib

.SUFFIXES : .o .f 

all:  ptor_driver.o 
	$(FC) $(FFLAGS) -o $(EXEC) ptor_driver.o  $(LLIB).a  $(EXTRA_LIBS)

.f.o : 
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -r *.o $(EXEC) $(LLIB).a
