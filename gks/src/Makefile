include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}
include ../install/make.ext.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG} ${GKSFLAGS}
else
   FFLAGS=${FOPT} ${GKSFLAGS}
endif

EXEC=gks_analysis
LLIB=gks_analysis_lib

OBJECTS = gks_var.o  \
          gks_analysis.o \
          gks_data_map.o   \
          gks_expprofiles.o \
          gks_tglf_map.o \
          gks_run.o     \
          gks_write_gksout.o       

.SUFFIXES : .o .f .f90

all: $(LLIB).a  gks_driver.o 
	$(FC) $(FFLAGS) -o $(EXEC) gks_driver.o  $(LLIB).a   $(EXTRA_LIBS)

gks_analysis_lib.a: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f.o :
	$(F77) $(FFLAGS) -I $(NETCDF_INC)  -c $<

clean:
	rm -f *.o *.lst $(LLIB).a $(EXEC)
	rm -f ${GACODE_ROOT}/modules/gks*.mod *mod
