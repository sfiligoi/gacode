include ${GYRO_DIR}/install/make.inc.${GYRO_PLAT}

EXEC = poincare

OBJ = input_data.o \
      runge_kutta.o \
      func.o \
      read_input.o \
      setup.o


.SUFFIXES : .o .f90 

all: $(OBJ) $(EXEC).o
	$(FC) -o $(EXEC) $(EXEC).o $(OBJ)

.f90.o :
	$(FC) $(FMATH) $(FDEBUG) -c $<

clean:  
	rm -f *.o $(EXEC) *~
	cd ${GYRO_DIR}/modules ; rm -f input_data.mod
