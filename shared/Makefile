include ../platform/build/make.inc.${GACODE_PLATFORM}

all:
	cd math ; make "OPT=${OPT}"
	cd landau ; make "OPT=${OPT}"
	cd harvest_client ; make
	cd nclass ; make "OPT=${OPT}"
	cd UMFPACK ; make "OPT=${OPT}"

clean:
	cd math ; make clean
	cd landau ; make clean
	cd harvest_client ; make clean
	cd nclass ; make clean
	cd UMFPACK ; make clean
