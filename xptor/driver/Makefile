include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

export EXTRA_LIBS = ${LMATH}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC= xptor
LLIB=xptor_lib

.SUFFIXES : .o .f 

all:  ptor_driver.o 
	$(FC) $(FFLAGS) -o $(EXEC) ptor_driver.o  $(LLIB).a  ${GACODE_ROOT}/tglf/src/tglf_lib.a $(EXTRA_LIBS) $(NETCDF) 

.f.o : 
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -r *.o $(EXEC) $(LLIB).a
