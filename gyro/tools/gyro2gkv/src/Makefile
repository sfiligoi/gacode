include ${GYRO_DIR}/install/make.inc.${GYRO_PLAT}

EXEC = gyro2gkv_exec

OBJECTS = gyro2gkv_globals.o \
          read_all.o \
          $(EXEC).o 

.SUFFIXES : .o .f90 

exec: $(OBJECTS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS) $(NETCDF)

.f90.o :
	$(FC) -I$(NETCDF_INC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -rf *.o *~ $(EXEC)
	cd ${GYRO_DIR}/modules ; rm -f gyro2gkv_globals.mod
