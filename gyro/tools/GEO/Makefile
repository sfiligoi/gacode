include ../../install/make.inc.${GYRO_PLAT}

EXEC=main
LLIB=../../GEO/GEO_lib.a

OBJECTS = $(EXEC).o

.SUFFIXES : .o .f90 

all: exec

exec: $(OBJECTS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS) $(LLIB)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o main
