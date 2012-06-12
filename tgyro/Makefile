#----------------------------
# Toplevel makefile for tgyro
#----------------------------

include ../shared/install/make.inc.${GACODE_PLATFORM}

export EXTRA_LIBS = \
       ${GACODE_ROOT}/neo/src/neo_lib.a \
       ${GACODE_ROOT}/neo/nclass/nclass_lib.a \
       ${GACODE_ROOT}/tglf/src/tglf_lib.a \
       ${GACODE_ROOT}/gyro/src/gyro_lib.a \
       ${GACODE_ROOT}/gyro/BLEND/BLEND_lib.a \
       ${GACODE_ROOT}/shared/EXPRO/EXPRO_lib.a \
       ${GACODE_ROOT}/shared/math/math_lib.a \
       ${GACODE_ROOT}/gyro/TRANSPOSE/TRANSP_lib.a \
       ${GACODE_ROOT}/gyro/SSUB/SSUB_lib.a \
       ${GACODE_ROOT}/shared/GEO/GEO_lib.a \
       ${GACODE_ROOT}/shared/UMFPACK/UMFPACK_lib.a \
       ${GACODE_ROOT}/gyro/GKEIGEN/GKEIGEN_lib.a \
       ${GACODE_ROOT}/tgyro/math/math_lib.a \
       ${GACODE_ROOT}/tgyro/qfm/qfm_lib.a \
       ${GACODE_ROOT}/gyro/vshdf5/gyrovshdf5_lib.a \
       ${GACODE_ROOT}/shared/vshdf5/vshdf5_lib.a \
       ${HDF5_LIBS} ${LMATH} ${LIPM}

all:
	gacode_getversion > .VERSION
	cd ${GACODE_ROOT}/shared ; make
	cd ${GACODE_ROOT}/neo ; make
	cd ${GACODE_ROOT}/tglf ; make
	cd ${GACODE_ROOT}/gyro ; make
	cd math ; make
	cd qfm ; make
	cd src ; make

clean:
	cd math ; make clean
	cd qfm ; make clean
	cd src ; make clean
	rm -rf .VERSION

.IGNORE:
