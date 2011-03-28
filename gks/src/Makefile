include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
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

.SUFFIXES : .o .f

all: $(LLIB).a  gks_driver.o 
	$(FC) $(FFLAGS) -o $(EXEC) gks_driver.o  $(LLIB).a   $(EXTRA_LIBS)

gks_analysis_lib.a: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f.o :
	$(FC) $(FMATH) $(FFLAGS) -I $(NETCDF_INC)  -c $<

clean:
	rm -f *.o $(LLIB).a $(EXEC)
	rm -f ${GACODE_ROOT}/modules/gks*.mod 
