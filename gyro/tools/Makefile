#-----------------------------
# Toplevel makefile for tools
#-----------------------------

include ../install/make.inc.${GYRO_PLAT}

all:
	cd ${GYRO_DIR}/math ; $(MAKE) "FFLAGS = ${FOPT}"
	cd ${GYRO_DIR}/GEO ; $(MAKE) "FFLAGS = ${FOPT}"
	cd ${GYRO_DIR}/tools/gyro2gkv/src ; $(MAKE) "FFLAGS = ${FOPT}"
	cd ${GYRO_DIR}/tools/iterdb_read/src ; $(MAKE) "FFLAGS = ${FOPT}"
	cd ${GYRO_DIR}/tools/iterdb_read/gato ; $(MAKE) "FFLAGS = ${FOPT}"

debug:
	cd ${GYRO_DIR}/math ; $(MAKE) "FFLAGS = ${FDEBUG}"
	cd ${GYRO_DIR}/GEO ; $(MAKE) "FFLAGS = ${FDEBUG}"
	cd ${GYRO_DIR}/tools/gyro2gkv/src ; $(MAKE) "FFLAGS = ${FDEBUG}"
	cd ${GYRO_DIR}/tools/iterdb_read/src ; $(MAKE) "FFLAGS = ${FDEBUG}"
	cd ${GYRO_DIR}/tools/iterdb_read/gato ; $(MAKE) "FFLAGS = ${FDEBUG}"

viz:
	cd ${GYRO_DIR}/tools/vtkgyro/torus ; $(MAKE)
	cd ${GYRO_DIR}/tools/spectra/src ; $(MAKE)

clean:
	cd gyro2gkv/src ; make clean
	cd iterdb_read/src ; make clean
	cd iterdb_read/gato ; make clean
	cd vtkgyro/torus ; make clean
	cd spectra/src ; make clean

help:
	cat README.help


