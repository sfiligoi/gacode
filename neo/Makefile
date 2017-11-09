#---------------------------
# Toplevel makefile for neo
#---------------------------

include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

export EXTRA_LIBS = \
       ${GACODE_ROOT}/shared/EXPRO/EXPRO_lib.a \
       ${GACODE_ROOT}/shared/math/math_lib.a \
       ${GACODE_ROOT}/shared/GEO/GEO_lib.a \
       ${GACODE_ROOT}/shared/UMFPACK/UMFPACK_lib.a \
       ${GACODE_ROOT}/shared/nclass/nclass_lib.a

#--------------------------------------------------------------------
# SOLVER_OPT can be set to any one of the following:
#     -- umf   : UMFPACK
#     -- slu   : SuperLU serial version
#     -- slumt : SuperLU multi-threaded version
# If SOLVER_OPT is not set to anything, the default option is 'umf'
#--------------------------------------------------------------------
SOLVER_OPT_DEFAULT = umf
ifeq ($(SOLVER_OPT),slu)
  BUILD_SUPERLU += slu
  EXTRA_LIBS += ${GACODE_ROOT}/shared/SUPERLU/SuperLU/libsuperlu.a
else ifeq ($(SOLVER_OPT),slumt)
  BUILD_SUPERLU += slumt
  EXTRA_LIBS += ${GACODE_ROOT}/shared/SUPERLU/SuperLU_MT/libsuperlu_mt.a
else
  export SOLVER_OPT := ${SOLVER_OPT_DEFAULT}
endif

all:
	gacode_getversion > .VERSION
	cd ${GACODE_ROOT}/shared ; make
	cd ${GACODE_ROOT}/shared/SUPERLU ; make ${BUILD_SUPERLU}
	cd src ; make
	cd ${GACODE_ROOT}/neo/tools/pneo ; make

clean:
	cd src ; make clean
	cd ${GACODE_ROOT}/neo/tools/pneo ; make clean

.IGNORE:
