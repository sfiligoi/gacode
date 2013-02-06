include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}
include ../install/make.ext.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG} ${XPTORFLAGS}
else
   FFLAGS=${FOPT} ${XPTORFLAGS}
endif

EXEC=xptor
LLIB=../driver/xptor_lib

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
	$(F77) $(FFLAGS) -I $(NETCDF_INC)  -c $<

.f90.o :
	$(FC) $(FFLAGS) -I $(NETCDF_INC)  -c $<


clean:
	rm -f *.o *.lst
