include install/make.inc.${GACODE_PLATFORM}

all:
	cd math ; make "OPT=${OPT}"
	cd GEO ; make "OPT=${OPT}"
	cd UMFPACK ; make "OPT=${OPT}"
	cd EXPRO ; make "OPT=${OPT}"
	cd fluxfit ; make "OPT=${OPT}"
	cd profiles_gen ; make "OPT=${OPT}"
	cd nclass ; make "OPT=${OPT}"

clean:
	cd math ; make clean
	cd GEO ; make clean
	cd UMFPACK ; make clean
	cd EXPRO ; make clean
	cd fluxfit ; make clean
	cd profiles_gen ; make clean
	cd nclass ; make clean
	cd python ; rm -f *.pyc	
