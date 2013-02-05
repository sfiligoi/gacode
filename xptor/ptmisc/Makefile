include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}
include ../install/make.ext.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG} ${XPTORFLAGS}
else
   FFLAGS=${FOPT} ${XPTORFLAGS}
endif

LLIB=../driver/xptor_lib

OBJECTS = ave_fields.o \
	  average7_2d.o \
          calc_stats.o \
          delay.o \
          derivedexpprofiles.o \
          derivedmodlocal.o \
          derivedmodprofiles.o \
          filter_dv.o \
          gridsmooth.o \
          init.o \
          ptorupdate.o \
          setup.o \
          timeprofiles_dv.o \
	  tpedscale.o \
          transfer_exp.o \
          update_src.o \
          write_ascii.o \
          write_ufile.o \
          write_netcdf.o \
          write_glf.o

.SUFFIXES : .o .f .F


lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS) 

.f.o :
	$(FC) $(FFLAGS) -c $<

write_netcdf.o:	write_netcdf.f
	$(FC) $(FFLAGS)  -I $(NETCDF_INC) -c $<

.F.o :
	$(FC) -c $(FFLAGS) $(SPECIAL) $<

clean:
	rm -f *.o *.lst 