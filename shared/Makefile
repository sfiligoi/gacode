include ../platform/build/make.inc.${GACODE_PLATFORM}

all:
	cd math ; make "OPT=${OPT}"
	cd nclass ; make "OPT=${OPT}"
	cd UMFPACK ; make "OPT=${OPT}"
ifdef FANN_ROOT
	cd ${GACODE_ROOT}/../neural && make
endif

clean:
	cd math ; make clean
	cd nclass ; make clean
	cd UMFPACK ; make clean
ifdef FANN_ROOT
	cd ${GACODE_ROOT}/../neural && make clean
endif

distclean:
	cd math ; make clean
