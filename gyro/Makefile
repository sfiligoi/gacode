include ../shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif
 
EXEC=gyro

export EXTRA_LIBS = \
       ${GACODE_ROOT}/shared/EXPRO/EXPRO_lib.a \
       ${GACODE_ROOT}/shared/math/math_lib.a \
       ${GACODE_ROOT}/shared/GEO/GEO_lib.a \
       ${GACODE_ROOT}/shared/UMFPACK/UMFPACK_lib.a \
       ${GACODE_ROOT}/gyro/BLEND/BLEND_lib.a \
       ${GACODE_ROOT}/gyro/TRANSPOSE/TRANSP_lib.a \
       ${GACODE_ROOT}/gyro/SSUB/SSUB_lib.a \
       ${GACODE_ROOT}/gyro/GKEIGEN/GKEIGEN_lib.a \
       ${LMATH}

all:
	cd ${GACODE_ROOT}/shared ; make
	cd BLEND ; make 
	cd TRANSPOSE ; make 
	cd SSUB ; make 
	cd src ; make
	cd GKEIGEN ; make
	cd peek ; make

	cd src ; $(FC) $(FFLAGS) -o $(EXEC) $(EXEC).o gyro_lib.a $(EXTRA_LIBS)

clean:
	cd BLEND ; make clean
	cd TRANSPOSE ; make clean
	cd SSUB ; make clean
	cd src ; make clean
	cd GKEIGEN ; make clean
	cd peek ; make clean

	cd src ; rm -f $(EXEC)

.IGNORE:
