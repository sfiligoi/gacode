#---------------------------
# Toplevel makefile for neo
#---------------------------

include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

#--------------------------------------------------------------------
# SOLVER_OPT can be set to any one of the following:
#     -- umf   : UMFPACK
#     -- slu   : SuperLU serial version
#     -- slumt : SuperLU multi-threaded version
# If SOLVER_OPT is not set to anything, the default option is 'umf'
#--------------------------------------------------------------------
export SOLVER_OPT ?= umf
ifeq ($(SOLVER_OPT),)
    export SOLVER_OPT := umf
endif

export EXTRA_LIBS = \
       ${GACODE_ROOT}/shared/EXPRO/EXPRO_lib.a \
       ${GACODE_ROOT}/shared/math/math_lib.a \
       ${GACODE_ROOT}/shared/GEO/GEO_lib.a \
       ${GACODE_ROOT}/shared/UMFPACK/UMFPACK_lib.a \
       ${GACODE_ROOT}/shared/nclass/nclass_lib.a

export SUPERLU_ROOT := ${GACODE_ROOT}/shared/SUPERLU
ifeq ($(SOLVER_OPT),slu)
   EXTRA_LIBS += ${SUPERLU_ROOT}/SuperLU/libsuperlu.a
else ifeq ($(SOLVER_OPT),slumt)
   EXTRA_LIBS += ${SUPERLU_ROOT}/SuperLU_MT/libsuperlu_mt.a
endif

all:
	gacode_getversion > .VERSION
	cd ${GACODE_ROOT}/shared ; make
	cd src ; make
	cd ${GACODE_ROOT}/neo/tools/pneo ; make

clean:
	cd src ; make clean
	cd ${GACODE_ROOT}/neo/tools/pneo ; make clean

.IGNORE:
