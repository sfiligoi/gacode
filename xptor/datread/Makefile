include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC=xptor
LLIB=../libs/xptor_lib

OBJECTS = readtestprofiles.o \
          readiterdb.o \
          readufiles.o \
          read_er.o \
          tenspline.o \
          gridsetup.o \
          datexch.o \
          datneo.o \
          datnclass.o \
          datmap.o \
          datavg.o \
          pressure.o \
          neutxp.o \
	  cdfread.o \
          u0read.o \
          u1read.o \
          u1tread.o \
          u2read.o \
          u2tread.o \
          u3read.o \
          cub_spline.o \
          average7_1d.o

.SUFFIXES : .o .f .f90


lib : $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS) 

.f.o :
	$(FC) $(FMATH) $(FFLAGS)  -c $<
.f90.o :
	$(FC) $(FMATH) $(FFLAGS)  -c $<


clean:
	rm -f *.o
