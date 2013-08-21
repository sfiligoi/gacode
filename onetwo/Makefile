#USE THIS MAKEFILE TO BUILD ONETWO On LOCAL ARCHITECTURES
# FORT=pgf90,FORT=ifort  are  valid,also OPT=optimize (default) or OPT=debug
#or OPT=debug  may be set on all machines
#gmake FORT=pgf95 OPT=optimize

#nbi_dimensions_mod,nbi_com_mod.f90,nbi_types.f90,ezspline_mod.f90,
#ezcdf.f90,ezcdf_genput.f90,ezcdf_opncls.f90,....,netcdf.inc,others ?? 
# are  picked up from nubeam area . 


 

SHELL     = /bin/sh
VERSION   = v5.4
GRID      = 129_201
NAME      = onetwo_$(GRID)
HOST      = $(shell uname -s)
ONETWO    = $(shell pwd)
SHARED_MODULES = $(ONETWO)/shared_modules
VPATH = $(SHARED_MODULES)
INSTALLDIR = ./ 
NUBEAM_12_BUILD = FALSE
UBUILD = FALSE
INCPATH = ./
NETCDF_INC =


ifeq ($(HOST),Linux)
        PORTLIB =	portlib
	IMSL = imsl_12
	CPPLF95LIB = 
	PGFLIBS   =
	ifeq ($(shell uname -n),isis1)
           include $(INCPATH)make.isis
	endif

	ifeq ($(shell uname -n),delphi2.gat.com)
           include $(INCPATH)make.delphi2
	endif
	ifeq ($(shell uname -n),lohan1.gat.com)
           #include $(INCPATH)make.lohans
           include $(INCPATH)make.lohan1
        endif
	ifeq ($(shell uname -n),lohan3.gat.com)
           #include $(INCPATH)make.lohans
           include $(INCPATH)make.lohan3
	endif
	ifeq ($(shell uname -n),benten.gat.com)
           #include $(INCPATH)make.lohans
           include $(INCPATH)make.lohan3
	endif
        ifeq ($(shell uname -n),hestia.gat.com)
           include $(INCPATH)make.zeus
        endif
	ifeq ($(shell uname -n),lohan4.gat.com)
           include $(INCPATH)make.lohan4
	endif
	ifeq ($(shell uname -n),lohan5.gat.com)
           include $(INCPATH)make.lohan6          
	endif

	ifeq ($(shell uname -n),head)
           include $(INCPATH)make.lohan6          
	endif

	ifeq ($(shell uname -n),node01.cluster)
           include $(INCPATH)make.lohan6          
	endif

	ifeq ($(shell uname -n),node06.cluster)
           include $(INCPATH)make.lohan6          
	endif

	ifeq ($(shell uname -n),zeus.gat.com)
           include $(INCPATH)make.zeus
	endif
	#ifeq ($(shell uname -n),hydra)
         #  include $(INCPATH)make.nfrc
	#endif
	ifeq ($(shell uname -n),jaguar1)
           include $(INCPATH)make.jaguar
	endif
	ifeq ($(shell uname -n),jaguar2)
           include $(INCPATH)make.jaguar
	endif
	ifeq ($(shell uname -n),jaguar3)
           include $(INCPATH)make.jaguar
	endif
	ifeq ($(shell uname -n),jaguar4)
           include $(INCPATH)make.jaguar
	endif
	ifeq ($(shell uname -n),jaguar5)
           include $(INCPATH)make.jaguar
	endif
	ifeq ($(shell uname -n),jaguar6)
           include $(INCPATH)make.jaguar
	endif
	ifeq ($(shell uname -n),aetna.gat.com)
           include $(INCPATH)make.aetna
	endif
	ifdef MACHINEOK
		#OPT can be changed on comand line:
		OPT=optimize
		#MPI = 1        # comment out this line if #MPI is not used
		#OPEN_MP = 1    #comment out this line if Open Mp is not used
                   #note NUBEAM modules are not compiled with openMp
                   #at this time (see compile_nubeam below,lines
                   #that start with C$ are activated with open mp directive)
		#PROFILE = 1     #comment out this line if 
				 #profileing  is not used
		#AUTO_PARA = 1   #comment out this line if auto 
				 #parallelization is not desired
		INCPATH	  = ./
	 	GCCLIB =  /usr/lib/libstdc++-3-libc6.2-2-2.10.0.so
                LISTING    = 
		ifeq ($(OPT),debug)
			OPTION = $(DEBUG) 
		else
			OPTION = $(OPTIMIZE) 
		endif
		ifndef OPEN_MP
			OMP =
		else
			#   do not include /usr/local/lib/libpgc.a for some reason or other
			#   it only allows us to use 1 cpu 
			OMP = -mp
			FTN90LIB =  \
				$(PGFLIBS)/lib/libpgftnrtl.a  \
				$(PGFLIBS)/lib/libpgthread.a \
				$(PGFLIBS)/lib/libpgmp.a
		endif

		ifndef PROFILE
			PROF = 
		else
			PROF = -Mprof=func
			PROF = -Mprof=lines
			FTN90LIB =  $(PGFLIBS)/lib/libpgc.a \
					$(PGFLIBS)/lib/libpgftnrtl.a \
					$(PGFLIBS)/lib/libpgthread.a \
					$(PGFLIBS)/lib/libpgmp.a \
					$(PGFLIBS)/lib/libpgnod_prof.a
              
		endif

		ifndef AUTO_PARA
			AUTP =
		else
			AUTP = -Mconcur -Minfo
			FTN90LIB = /usr/local/lib/libpgftnrtl.a \
				 /usr/local/lib/libpgthread.a \
				  /usr/local/lib/libpgmp.a
		endif


		LOAD      = $(FC) $(LINK)  -o $(NAME) # -m
		ARCHIVE   = ar rv
	endif

endif   # over HOSTS 


PROTECT   = chmod 755
DELETE    = rm -f
RENAME    = mv
include $(INCPATH)/$(MAKE1)
include $(INCPATH)/make.sources
include $(INCPATH)/$(MAKE2)
 
#note to me portlib was changed !!!!!!!!!!!!!!!!
#NOTE***********************************************************************
#
#Now, ONETWO does not call the NUBEAM related libraries anymore.
#
#Note carefully that f90 modules, say a.f90, when compiled produce a.mod
#and a.o . You must link with a.o, but a.mod is also required. 
#failure to link with a.o does NOT produce a link error!!
#but the module info will not be there as expected - Linux with pgf90 
#the search order of the libraries given below does not agree
#with the order posted in nubeam readme but is necessary to get linkup!!!!!
#****************************************************************************


ifeq ($(strip $(FORT)),ifort)
##LIBRARIES = $(LIBS)/$(PORTLIB).a $(LIBS)/cortlib.a  \
##			$(LIBS)/$(IMSL).a \
##		        $(LIBS)/method_lines_library.a  $(LAPACK) \
##			$(UFILE_LIBRARIES) $(NETCDF_LIB)  $(TSPLINE_LIB)/$(SPLINE_LIB)
LIBRARIES = $(LIBS)/$(PORTLIB).a $(LIBS)/cortlib.a \
		        $(LIBS)/method_lines_library.a  $(LAPACK) \
			$(UFILE_LIBRARIES) $(NETCDF_LIB)  $(TSPLINE_LIB)/$(SPLINE_LIB)
endif
ifeq ($(strip $(FORT)),pgf95)
##LIBRARIES = $(LIBS)/$(PORTLIB).a $(LIBS)/cortlib.a   $(NETCDF_LIB) \
##			$(LIBS)/$(IMSL).a \
##		        $(LIBS)/method_lines_library.a  $(LAPACK) \
##			$(UFILE_LIBRARIES)  $(TSPLINE_LIB)/$(SPLINE_LIB)
LIBRARIES = $(LIBS)/$(PORTLIB).a $(LIBS)/cortlib.a   $(NETCDF_LIB) \
		        $(LIBS)/method_lines_library.a  $(LAPACK) \
			$(UFILE_LIBRARIES)  $(TSPLINE_LIB)/$(SPLINE_LIB)
endif
.SUFFIXES:


all: $(NAME) data_files
	@echo
	@echo ">>  `date '+%a %d-%h-%y %r'`  `pwd`  `uname -mns`  $(LOGNAME)"
	@echo
	@echo make of $(NAME) in `pwd` on `hostname` completed
#	@echo $(UFILE_LIBRARIES)
#	@echo $(NETCDF_LIB)

$(NAME):           $(OBJECTS) $(OBJECTS_P) $(LIBRARIES)
	@echo "loading .... "
	$(LOAD)    $(OBJECTS)  $(LIBRARIES)    $(FTN90LIB)
	$(PROTECT) $(NAME)
#	cp $(NAME) $(INSTALLDIR)

#JMP START

ifdef PARALLEL_12 #------------------------------------------------------------

mpi12.o: mpi12.F90 
	$(COMPILE) mpi12.F90

onetwo_parallel.o:   onetwo_parallel.F90 mpi12.o param.o
	$(COMPILE) onetwo_parallel.F90

onetwo_slave.o:    onetwo_slave.F90 mpi12.o param.o
	$(COMPILE) onetwo_slave.F90

glf_slave.o:    glf_slave.F90 mpi12.o
	$(COMPILE) glf_slave.F90

endif #-------------------------------------------------------------------------

cray001.o:         cray001.F   $(INCLUDES_cray001)
	$(COMPILE) cray001.F

cray101.o:         prl_mod.o cray341.o cray101.f   $(INCLUDES_cray101)
	$(COMPILE) cray101.f

cray102.o:         cray102.f   $(INCLUDES_cray102) aid_newton.o 
	$(COMPILE) cray102.f  

cray202.o:         cray202.f   $(INCLUDES_cray202)
	$(COMPILE) cray202.f

cray203.o:         cray203.f   $(INCLUDES_cray203)
	$(COMPILE) cray203.f

cray204.o:         cray204.f   $(INCLUDES_cray204)
	$(COMPILE) cray204.f

cray205.o:         cray205.f   $(INCLUDES_cray205)
	$(COMPILE) -I$(NETCDF_INC)  cray205.f

cray206.o:         cray206.f   $(INCLUDES_cray206)
	$(COMPILE)  -I$(NETCDF_INC) cray206.f

cray207.o:         cray207.f   $(INCLUDES_cray207)
	$(COMPILE) cray207.f

cray208.o:         cray208.f   $(INCLUDES_cray208)
	$(COMPILE) cray208.f

cray209.o:         cray209.f   $(INCLUDES_cray209)
	$(COMPILE) cray209.f

cray302.o:         cray302.f   $(INCLUDES_cray302)
	$(COMPILE) cray302.f

cray303.o:         cray303.f   $(INCLUDES_cray303)
	$(COMPILE) cray303.f

cray304.o:         cray304.f   $(INCLUDES_cray304)
	$(COMPILE) cray304.f

cray306.o:         cray306.f cray341.o    $(INCLUDES_cray306)
	$(COMPILE) cray306.f

cray307.o:         cray307.F   $(INCLUDES_cray307) $(INCLUDES_P_cray307)
	$(COMPILE) cray307.F

cray308.o:         cray308.f   $(INCLUDES_cray308)
	$(COMPILE) cray308.f

cray309_nubeam.o:         cray309_nubeam.f   $(INCLUDES_cray309)
	$(COMPILE) cray309_nubeam.f

cray311.o:         cray311.f   $(INCLUDES_cray311)
	$(COMPILE) cray311.f

cray321u.o:         cray321u.f   $(INCLUDES_cray321)
	$(COMPILE) cray321u.f

cray331.o:         cray331.f   $(INCLUDES_cray331)
	$(COMPILE) -I$(NETCDF_INC)  cray331.f

cray341.mod:	cray341.o
cray341.o:         prl_mod.o spec_kind_mod.o cray341.f90   $(INCLUDES_cray341)
	$(COMPILE) cray341.f90

cray352.o:         cray352.f   $(INCLUDES_cray352)
	$(COMPILE) cray352.f

cray401.o:         cray401.f   $(INCLUDES_cray401)
	$(COMPILE) cray401.f

cray402.o:         cray402.f   $(INCLUDES_cray402)
	$(COMPILE) cray402.f

cray403.o:         cray403.f   $(INCLUDES_cray403)
	$(COMPILE) cray403.f

cray404.o:         cray404.f   $(INCLUDES_cray404)
	$(COMPILE) cray404.f

cray501.o:         cray501.f   $(INCLUDES_cray501)
	$(COMPILE) cray501.f

callglf2d.o:        callglf2d.f   glf.m
	$(COMPILE) callglf2d.f

writeglf2d.o:      writeglf2d.f  
	$(COMPILE) writeglf2d.f

glf2d.o:           glf2d.F   $(INCLUDES_glf2d)
	$(COMPILE) glf2d.F

r8tomsqz.o:        r8tomsqz.F f77_dcomplx.h
	$(COMPILE)  r8tomsqz.F

cray504.o:         cray504.f   $(INCLUDES_cray504)
	$(COMPILE) cray504.f

cray601.o:         cray601.f   $(INCLUDES_cray601) $(INCLUDES_MODULES)
	$(COMPILE) cray601.f

cray700.o:         cray700.f   $(INCLUDES_cray700)  $(INCLUDES_MODULES)
	#echo $(INCLUDES_cray700) $(INCLUDES_MODULES)
	$(COMPILE) cray700.f

allocate_plasma_freq.o:	allocate_plasma_freq.f90  $(INCLUDES_pl_freq)
	$(COMPILE)  -DONETWO -I$(NETCDF_INC)  $< -o  $@ 

cfactrs.mod: confinement_factors.o
confinement_factors.o:	confinement_factors.f90 $(INCLUDES_cfactors)
	$(COMPILE) confinement_factors.f90

check_en_bc.o:	check_en_bc.f90 nrtype.o param.o numbrs.o bd_condtn.mod
	$(COMPILE) check_en_bc.f90

bicubic_spline.o:         bicubic_spline.f90  nrtype.o io_gcnmp.o error_handler.o
	$(COMPILE)  $(GFLAGS)  $< 



rw_iterdb.o:	rw_iterdb.F90 $(INCLUDES_iterdb_gcnmp)
	$(COMPILE)  -DONETWO -I$(NETCDF_INC)  $< -o  $@ 

#imfit_iterdb.o:	imfit_iterdb.f90 $(INCLUDES_IMFIT_ITERDB)
#	$(COMPILE) imfit_iterdb.f90

fdyimp.mod:	fdyimp.o
fdyimp.o:	fdyimp.f90 $(INCLUDES_FDYIMP)
	$(COMPILE) $< -o $@ 

fusion_gcnmp.o: fusion_gcnmp.f90  nrtype.o  error_handler.o consts.o grid_class.o \
                  ions_gcnmp.o dep_var_gcnmp.o source_terms_gcnmp.o
	$(COMPILE) $< -o $@ 

gcnmp_version.mod:	gcnmp_version.o
gcnmp_version.o:	gcnmp_version.f90 $(INCLUDES_GCNMP_VERSION)
	$(COMPILE) $< -o $@ 

gcnmp_input.mod:	gcnmp_input.o param.o
gcnmp_input.o:	gcnmp_input.f90 nrtype.o solcon_gcnmp.o  set_12_gcnmp_vars.o
	$(COMPILE) gcnmp_input.f90

gcnmp_interface.mod:	gcnmp_interface.o 
gcnmp_interface.o : gcnmp_interface.f90 gcnmp_namelist.o gcnmp_input.o
	$(COMPILE) gcnmp_interface.f90

gcnmp_namelist.o : gcnmp_namelist.F90 $(INCLUDES_netcdf_read) \
		  bc_values_gcnmp.o fusion_gcnmp.o tglfin.o \
		  tglf_input_defaults.o turb_tport.o mmm_in.o \
		  forcebal_data_mod.o
	$(COMPILE) $< -o $@ 

gcnmp_namelist_defaults.o : gcnmp_namelist_defaults.f90  \
		  bc_values_gcnmp.o nrtype.o source_terms_gcnmp.o neo_tport.o
	$(COMPILE) $< -o $@ 

gpsi.mod:	gpsi.o
gpsi.o : gpsi.f90 $(INCLUDES_netcdf_read) bc_values_gcnmp.o
	$(COMPILE) gpsi.f90

island.mod:	island.o
island.o : island.f90 nrtype.o
	$(COMPILE) island.f90

rw_iterdb_netcdf.o : rw_iterdb_netcdf.F90 $(INCLUDES_netcdf_read)
	$(COMPILE) -I$(NETCDF_INC) $< -o $@

rw_iterdb_netcdf_mod.o : rw_iterdb_netcdf_mod.F90 $(INCLUDES_netcdf_read)
	$(COMPILE) -I$(NETCDF_INC) $< -o $@

iterdb_v1.29.o:	iterdb_v1.29.f  $(INCLUDES_iterdb_v129.f)
	$(COMPILE) iterdb_v1.29.f 

method_lines.o:    method_lines.f $(INCLUDES_method_lines)
	$(COMPILE) method_lines.f

ohacd.o:         ohacd.f90   $(INCLUDES_ohacd)
	$(COMPILE) ohacd.f90

load_12.o:         load_12.f90   transp.o
	$(COMPILE) load_12.f90

load_restart_12.o:         load_restart_12.f90    restore_12.o tordlrot.mod
	$(COMPILE) load_restart_12.f90

kapisn.o:	kapisn.F90 $(INCLUDES_kapisn)
	$(COMPILE) $< -o  $@ 

restore_12.mod:		restore_12.o 
restore_12.o:           restore_12.f90 nrtype.o  
	$(COMPILE) restore_12.f90

monitr.mod:		monitr.o  param.o
monitr.o:           monitr.f90
	$(COMPILE) monitr.f90


nrtype.mod:		nrtype.o 
nrtype.o:           nrtype.F90
	@echo  FORT=$(FORT)     
	@echo COMPILE   = $(COMPILE)
	@echo USEMPI    = $(USEMPI)
	@echo MPI2_EXT  = $(MPI2_EXT)
	@echo MPI_INC   = $(MPI_INC)
	@echo ">>  `date '+%a %d-%h-%y %r'`  `pwd`  `uname -mns`  $(LOGNAME)"
	$(COMPILE) $<  -o $@

rad_loss.mod:	rad_loss.o nrtype.o 
rad_loss.o:           rad_loss.f90
	$(COMPILE) rad_loss.f90

#run_gcnmp.mod:	run_gcnmp.o 
#run_gcnmp.o:           run_gcnmp.f90 nrtype.o gcnmp_interface.o io.o
#	$(COMPILE) run_gcnmp.f90

shot_info.mod:	shot_info.o
shot_info.o:	shot_info.f90 nrtype.o
	$(COMPILE) $< -o $@ 


grid_class.mod:	grid_class.o 
grid_class.o:	grid_class.F90 nrtype.o  tspline.o  error_handler.o vector_class.o nrtype.o solcon_gcnmp.o
	$(COMPILE) $< -o $@ 

string_util.mod:	string_util.o  
string_util.o:    string_util.F90 mesh.o param.o
	$(COMPILE)  -DONETWO  $< -o  $@
 
replace_imsl.mod:	replace_imsl.o 
replace_imsl.o:		replace_imsl.F90 nrtype.o m_refsor.o \
			bicube.o cubic_spline_class.o  \
			bicubic_spline.o
	$(COMPILE)   $< -o  $@ 
#	$(COMPILE)  -DTESTING   $< -o  $

m_refsor.mod:	m_refsor.o
m_refsor.o:	m_refsor.f90 nrtype.o 
	$(COMPILE)   $< -o  $@ 


thermalization.o:	thermalization.f90 nrtype.o consts.o constnts.o \
			error_handler.o io_gcnmp.o nub2.o
	$(COMPILE) $< -o $@ 

aid_newton.mod:         aid_newton.o  
aid_newton.o:         aid_newton.f90   
	$(COMPILE) aid_newton.f90



bc_values_gcnmp.mod:         bc_values_gcnmp.o 
bc_values_gcnmp.o:         bc_values_gcnmp.F90   nrtype.o solcon_gcnmp.o \
		    	   getargs.o get_next_io_unit.o io_gcnmp.o grid_class.o
	$(COMPILE) -DONETWO  $< -o $@ 

cer.mod:         cer.o param.o nrtype.mod  mhdpar.mod
cer.o:         cer.f90  param.o nrtype.o mhdpar.o  
	$(COMPILE) cer.f90


consts.mod:         consts.o 
common_constants.mod:   consts.o
consts.o:         consts.f90   nrtype.o 
	$(COMPILE) $< -o $@

dep_var_gcnmp.mod:         dep_var_gcnmp.o 
dep_var_gcnmp.o:         dep_var_gcnmp.f90   nrtype.o common_constants.mod
	$(COMPILE) $< -o $@ 

etc.mod:         etc.o  
etc.o:         etc.f90   param.o  
	$(COMPILE) etc.f90


iterdbmd.mod:         iterdbmd.o   nrtype.mod
iterdbmd.o:         iterdbmd.f90   nrtype.o
	$(COMPILE) iterdbmd.f90


iterdbmd_gcnmp.mod:         iterdbmd_gcnmp.o   nrtype.mod
iterdbmd_gcnmp.o:         iterdbmd_gcnmp.f90   nrtype.o grid_class.o solcon_gcnmp.o
	$(COMPILE) $< -o $@ 


mixcom.mod:         mixcom.o 
mixcom.o:         mixcom.f90   
	$(COMPILE) mixcom.f90

fast_ion_diffusion.mod:    fast_ion_diffusion.o 
fast_ion_diffusion.o:        fast_ion_diffusion.f90   
	$(COMPILE) fast_ion_diffusion.f90


psig.mod:         psig.o 
psig.o:         psig.f90  param.o mhdpar.o
	$(COMPILE) psig.f90


terminate.mod:         terminate.o 
terminate.o:         terminate.f90   
	$(COMPILE) terminate.f90


turb_tport.o :	turb_tport.f90  nrtype.o solcon_gcnmp.o 
	$(COMPILE) $< -o $@ 
weiland.mod:	weiland.o
weiland.o : weiland.f90 nrtype.o
	$(COMPILE) weiland.f90

zeffcom.mod:	zeffcom.o
zeffcom.o : zeffcom.f90 nrtype.o
	$(COMPILE) zeffcom.f90


allocate_err.mod:         allocate_err.o  
allocate_err.o:         allocate_err.f90   
	$(COMPILE)   $< -o  $@



com.mod:         com.o   
com.o:         com.f90   
	$(COMPILE) com.f90


compact_schemes.o:	compact_schemes.f90 nrtype.o error_handler.o io_gcnmp.o
	$(COMPILE) $< -o $@ 

copy.mod:         copy.o    
copy.o:         copy.f90    
	$(COMPILE) copy.f90


colrate.mod:         colrate.o    
colrate.o:         colrate.f90    
	$(COMPILE) colrate.f90


file_proc.mod:	file_proc.o 
file_proc.o:	file_proc.f90 nrtype.o solcon.o io.mod echdat_module.mod
	$(COMPILE) $< -o $@ 

flx.mod:         flx.o   
flx.o:         flx.f90   
	$(COMPILE) flx.f90



glf23.mod:         glf23.o  
glf23.o:         glf23.f90   
	$(COMPILE) glf23.f90


glf23_gcnmp.mod:         glf23_gcnmp.o  
glf23_gcnmp.o:         glf23_gcnmp.f90   
	$(COMPILE) $< -o $@ 

glf_tport_gcnmp.mod:         glf_tport_gcnmp.o  
glf_tport_gcnmp.o:         glf_tport_gcnmp.F90 tspline.o glf23_gcnmp.o \
                        neutral_beams.o curden_terms.o 
	$(COMPILE) $< -o $@


mesh.mod:         mesh.o  
mesh.o:         mesh.f90   param.o
	$(COMPILE) mesh.f90

neutral_beams.o:	neutral_beams.F90 nf_param.o	nrtype.o ions_gcnmp.o \
			plasma_prop.o nrtype.o consts.o grid_class.o nf_param.o
	$(COMPILE) $< -o $@


neo_tport.mod:         neo_tport.o   
neo_tport.o:         neo_tport.f90  neutral_beams.o source_terms_gcnmp.o 
	$(COMPILE) $< -o $@

spec_kind_mod.mod:       spec_kind_mod.o
spec_kind_mod.o:         spec_kind_mod.f
	$(COMPILE) $< -o $@

prl_mod.mod:          prl_mod.o
prl_mod.o:            prl_mod.f90 spec_kind_mod.o 
	$(COMPILE) -I$(NETCDF_INC) $< -o $@

ode_pkg.mod:         ode_pkg.o  
ode_pkg.o:         ode_pkg.f90   
	$(COMPILE) ode_pkg.f90

#ifs.mod:         ifs.o  
ifs.o ifs.mod:         ifs.f90   mhdpar.o
	$(COMPILE) ifs.f90


Statefile_proc.o:         Statefile_proc.f90 numbrs.o nrtype.o \
			  iterdbmd.o iterdbmd_gcnmp.o  gcnmp_input.o
	$(COMPILE) Statefile_proc.f90


sourc.mod:         sourc.o  
sourc.o:         sourc.f90   
	$(COMPILE) sourc.f90


verbose.mod:         verbose.o  
verbose.o:         verbose.f90   
	$(COMPILE) verbose.f90



adaptive.mod:         adaptive.o  
adaptive.o:         adaptive.f90   param.o
	$(COMPILE) adaptive.f90



eq_pfprim.mod:          eq_pfprim.o
eq_pfprim.o:		eq_pfprim.f90  io.mod machin.mod
	$(COMPILE) eq_pfprim.f90


neut.mod:         neut.o
neut.o:		neut.f90  
	$(COMPILE) neut.f90


precision_mod.mod:         precision_mod.o
precision_mod.o:		precision_mod.f90  
	$(COMPILE) precision_mod.f90


neucg.mod:         neucg.o
neucg.o:		neucg.f90  
	$(COMPILE) neucg.f90


paleocl.mod:            paleocl.o 
paleocl.o:		paleocl.f90  flx.o
	$(COMPILE) paleocl.f90

pelcom.mod:            pelcom.o 
pelcom.o:	pelcom.f90  nrtype.o
	$(COMPILE) pelcom.f90

plasma_prop.o:	plasma_prop.f90  vector_class.o ions_gcnmp.o
	$(COMPILE) $< -o $@ 

kinetic_efit.mod:         kin_efit.o  
kin_efit.o:		kin_efit.f90  $(INCLUDES_kin_efit)
	$(COMPILE) kin_efit.f90


staebler.mod:         staebler.o 
staebler.o:         staebler.f90 param.o  
	$(COMPILE) staebler.f90


tcoef.mod:         tcoef.o param.o  nrtype.o
tcoef.o:         tcoef.F90   
	$(COMPILE)  -DONETWO    $<  -o $@


rhog.mod:         rhog.o param.o  
rhog.o:         rhog.f90   
	$(COMPILE) rhog.f90


mhdcom.mod:         mhdcom.o    mhdpar.mod
mhdcom.o:         mhdcom.f90   mhdpar.o
	$(COMPILE) mhdcom.f90

solcon.mod:         solcon.o  
solcon.o:         solcon.f90  param.o 
	$(COMPILE) solcon.f90


solcon_gcnmp.mod:         solcon_gcnmp.o 
solcon_gcnmp.o:         solcon_gcnmp.f90 io_gcnmp.o  ions_gcnmp.o \
			param.o consts.o error_handler.o dep_var_gcnmp.o
	$(COMPILE) $< -o $@ 


source_terms_gcnmp.mod:         source_terms_gcnmp.o param.o 
source_terms_gcnmp.o:         source_terms_gcnmp.f90 io_gcnmp.o  ions_gcnmp.o \
                                solcon_gcnmp.o grid_class.o
	$(COMPILE) $< -o $@ 
ifeq ($(HOST),HP-UX)
	cp SOURCE_TERMS_GCNMP.mod source_terms_gcnmp.mod
endif
soln.mod:         soln.o 
soln.o:         soln.f90 param.o
	$(COMPILE) soln.f90


soln2d.mod:         soln2d.o
soln2d.o:         soln2d.f90 mhdpar.o param.o
	$(COMPILE) soln2d.f90
ifeq ($(HOST),HP-UX)
	cp SOLN2D.mod soln2d.mod
endif

tglfin.mod:         tglfin.o 
tglfin.o:	tglfin.f90 tglfv1.94_modules.o   nrtype.o tglf_modules.o
	$(COMPILE) $< -o $@ 

tglf_modules.o:	tglf_modules.f90
	$(COMPILE) $< -o $@ 

tglfv1.94_modules.o:	tglfv1.94_modules.F90
	$(COMPILE)  -DONETWO    $<  -o $@

tglf_input_defaults.mod:         tglf_input_defaults.o  
tglf_input_defaults.o:	tglf_input_defaults.f90  nrtype.o tglfin.o solcon_gcnmp.mod
	$(COMPILE) $< -o $@ 

fusion.mod:         fusion.o param.o
fusion.o:	fusion.f90 param.mod
	$(COMPILE) fusion.f90
ifeq ($(HOST),HP-UX)
	cp FUSION.mod fusion.mod
endif

ions.mod:         ions.o param.o
ions.o:	ions.f90 param.mod
	$(COMPILE) ions.f90

ions_gcnmp.mod:         ions_gcnmp.o param.o
ions_gcnmp.o:	ions_gcnmp.f90 param.o
	$(COMPILE) $< -o $@  

paleo.o: paleo.f90
	$(COMPILE) $<
test_paleo.o: test_paleo.f90 paleo.o terminate.o
	$(COMPILE) $<
test_paleo: paleo.o test_paleo.o terminate.o
	${LOAD} -o $@ paleo.o test_paleo.o terminate.o


forcebal_data_mod.mod:  forcebal_data_mod.o spec_kind_mod.mod
forcebal_data_mod.o: forcebal_data_mod.f90 spec_kind_mod.o
	$(COMPILE) $<

machin.mod:         machin.o param.mod
machin.o:	machin.f90 param.mod
	$(COMPILE) machin.f90


mhdbcdtn.mod:         mhdbcdtn.o param.mod
mhdbcdtn.o:	mhdbcdtn.f90 param.mod
	$(COMPILE) mhdbcdtn.f90

nub.mod:         nub.o param.mod
nub.o:	nub.F90 param.o
	$(COMPILE) -DONETWO $< -o $@


nub2.mod:         nub2.o param.mod
nub2.o:	nub2.f90 param.o
	$(COMPILE) nub2.f90


nub3.mod:         nub3.o param.mod
nub3.o:	nub3.f90 param.o
	$(COMPILE) nub3.f90


nub4.mod:         nub4.o param.mod
nub4.o:	nub4.f90 param.mod
	$(COMPILE) nub4.f90


oloss_dat.mod:         oloss.o param.mod
oloss.o:	oloss.f90 param.mod
	$(COMPILE) oloss.f90

ename.mod:	ename.o       
ename.o:	ename.f90 solcon.o yoka.o
	$(COMPILE) ename.f90


m_lines.mod:          m_lines.o
m_lines.o:          m_lines.f90   
	$(COMPILE)  m_lines.f90


bicube.mod:          bicube.o
bicube.o:          bicube.f90   cubic_spline.mod
	$(COMPILE)  $< -o  $@


flxav.mod:          flxav.o
flxav.o:          flxav.f90   
	$(COMPILE)  flxav.f90


shapctr.mod:          shapctr.o
shapctr.o:          shapctr.f90  mhdpar.o 
	$(COMPILE)  shapctr.f90


metrics.mod:          metrics.o
metrics.o:          metrics.f90   mhdpar.o
	$(COMPILE)  metrics.f90


neo2d.mod:          neo2d.o
neo2d.o:          neo2d.f90   
	$(COMPILE)  neo2d.f90

 
neo2dp.mod:          neo2dp.o
neo2dp.o:          neo2dp.f90   mhdpar.o
	$(COMPILE)  neo2dp.f90


ext_prog_info.mod:         ext_prog_info.o  
ext_prog_info.o:	ext_prog_info.f90  mhdpar.o io.o
	$(COMPILE) ext_prog_info.f90


toq_12_interface.mod:	toq_12_interface.o 
toq_12_interface.o:	toq_12_interface.f90 ext_prog_info.mod
	$(COMPILE) toq_12_interface.f90


toq_12.o:	toq_12.f90 $(INCLUDES_toq_12)
	$(COMPILE) toq_12.f90

limiter.mod:         limiter.o
limiter.o:	     limiter.f90
	$(COMPILE)   limiter.f90 


mhdpar.mod:         mhdpar.o
mhdpar.o:	mhdpar.f90 param.o
	$(COMPILE) $< -o $@ 


contour.mod:         contour.o
contour.o:	contour.f90
	$(COMPILE)   contour.f90 

Collect_stats.mod:	error_handler.o io_gcnmp.o nrtype.o
Collect_stats.o:	Collect_stats.f90
	$(COMPILE) $< -o $@   


tdem.mod:         tdem.o 
tdem.o:	          tdem.f90 contour.o limiter.o mhdpar.mod limiter.mod contour.mod param.mod
	$(COMPILE)   tdem.f90


param.mod: param.o         
param.o:	param.f90
	$(COMPILE) $< -o $@ 

bd_condtn.mod: bd_condtn.o         
bd_condtn.o:	bd_condtn.f90 nrtype.o consts.o param.o
	$(COMPILE)   bd_condtn.f90

io.mod: io.o       
io.o:	io.f90  param.o
	$(COMPILE)   io.f90

io_gcnmp.mod: io_gcnmp.o       
io_gcnmp.o:	io_gcnmp.f90 
	$(COMPILE) $< -o $@  

flags.mod: flags.o       
flags.o:	flags.f90  param.o
	$(COMPILE)   flags.f90


fast_ion_data_gcnmp.mod: fast_ion_data_gcnmp.o       
fast_ion_data_gcnmp.o:	fast_ion_data_gcnmp.f90 nrtype.o
	$(COMPILE) $< -o $@

neutral_data.mod: neutral_data.o       
neutral_data.o:	neutral_data.f90 nrtype.o
	$(COMPILE) $< -o $@ 

error_handler.mod:	error_handler.o io_gcnmp.o 
error_handler.o:	error_handler.F90 MPI_data.o
	$(COMPILE) $< -o $@   

extra.mod: extra.o       
extra.o:	extra.f90  param.o
	$(COMPILE)   extra.f90
ifeq ($(HOST),HP-UX)
	cp EXTRA.mod extra.mod
endif

getargs.mod: getargs.o       
getargs.o:	getargs.f90  param.o gcnmp_version.o error_handler.o
	$(COMPILE)   getargs.f90
getargs.mod: getargs.o   

get_next_io_unit.mod: get_next_io_unit.o       
get_next_io_unit.o:	get_next_io_unit.F90  param.o error_handler.o
	$(COMPILE)     get_next_io_unit.F90

geom.mod:	geom.o
geom.o:	geom.f90  param.o
	$(COMPILE)   geom.f90
ifeq ($(HOST),HP-UX)
	cp GEOM.mod geom.mod
endif

MPI_data.mod:	MPI_data.o
MPI_data.o:	MPI_data.F90  nrtype.o consts.o
	$(COMPILE) $< -o $@   

yoka.mod: yoka.o       
yoka.o:	yoka.f90  param.o
	$(COMPILE)   yoka.f90
ifeq ($(HOST),HP-UX)
	cp YOKA.mod yoka.mod
endif

mcgo.mod: mcgo.o       
mcgo.o:	mcgo.f90  param.o
	$(COMPILE)   mcgo.f90
ifeq ($(HOST),HP-UX)
	cp MCGO.mod mcgo.mod
endif

mhdgrid.mod: mhdgrid.o       
mhdgrid.o:	mhdgrid.f90  param.o mhdpar.o
	$(COMPILE)   mhdgrid.f90
ifeq ($(HOST),HP-UX)
	cp MHDGRID.mod mhdgrid.mod
endif

rf.mod: rfmod.o param.o       
rfmod.o:	rfmod.f90 soln.o param.o ename.o solcon.o io.o ions.o \
	 	param.o fusion.o nub.o  nub2.o numbrs.o
	$(COMPILE)   rfmod.f90
ifeq ($(HOST),HP-UX)
	cp RF.mod rfmod.mod
endif

echdat_module.mod:	echdat_module.o param.o io.o rfmod.o   
echdat_module.o:	echdat_module.F90  param.o  io.o rfmod.o \
			 error_handler.o
	$(COMPILE) -DONETWO -I$(NETCDF_INC) $< -o $@ 

numbrs.mod: numbrs.o       
numbrs.o: param.mod
	$(COMPILE)   numbrs.f90
ifeq ($(HOST),HP-UX)
	cp NUMBRS.mod numbrs.mod
endif

nonlin.mod: nonlin.o       
nonlin.o:
	$(COMPILE)   nonlin.f90
ifeq ($(HOST),HP-UX)
	cp NONLIN.mod NONLIN.mod
endif

tfact.mod: tfact.o  
tfact.o:	tfact.f90
	$(COMPILE)   tfact.f90
ifeq ($(HOST),HP-UX)
	cp TFACT.mod tfact.mod
endif

xptor_sim.mod: xptor_sim.o  
xptor_sim.o:	xptor_sim.f90
	$(COMPILE)   xptor_sim.f90
ifeq ($(HOST),HP-UX)
	cp XPTOR_SIM.mod xptor_sim.mod
endif

d01bff_dat.mod: d01bff_dat.o  
d01bff_dat.o:	d01bff_dat.f90
	$(COMPILE)   d01bff_dat.f90
ifeq ($(HOST),HP-UX)
	cp D01BFF_DAT.mod d01bff_dat.mod
endif

cpub_dat.mod: cpub_dat.o  
cpub_dat.o:	cpub_dat.f90
	$(COMPILE)   cpub_dat.f90
ifeq ($(HOST),HP-UX)
	cp CPUB_DAT.mod cpub_dat.mod
endif
tordlrot.mod:         tordlrot.o  
tordlrot.o:           tordlrot.f90   param.o
	$(COMPILE) tordlrot.f90
ifeq ($(HOST),HP-UX)
	cp TORDLROT.mod tordlrot.mod
endif

tmpcom.mod:         tmpcom.o  
tmpcom.o:           tmpcom.f90   param.o
	$(COMPILE) tmpcom.f90
ifeq ($(HOST),HP-UX)
	cp TMPCOM.mod tmpcom.mod
endif

cubic_spline.mod: cubic_spline_class.o
cubic_spline_class.o:	cubic_spline_class.f90 nrtype.o common_constants.mod
	$(COMPILE) $<

constnts.mod:         constnts.o  
constnts.o:           constnts.f90 
	$(COMPILE) constnts.f90
 
pnpoly.o:	pnpoly.f90 
	$(COMPILE)  pnpoly.f90

smthspline.mod:         nrtype.mod smthspline.o
smthspline.o:           nrtype.o smthspline.f90 
	$(COMPILE) smthspline.f90

mmm95_ntcc.o:           mmm95_ntcc.f 
	$(COMPILE) mmm95_ntcc.f

mmm95_mod.mod:         mmm95_mod.o
mmm95_mod.o:           mmm95_mod.f90 
	$(COMPILE) mmm95_mod.f90 

set_12_gcnmp_vars.o:	set_12_gcnmp_vars.f90 tglfin.o $(INCLUDES_netcdf_read) zen_state.o \
			fusion_gcnmp.o plasma_prop.o P_Nfreya_12_interface.o shot_info.o \
	                glf23_gcnmp.o yoka.o flags.o tcoef.o gpsi.mod allocate_plasma_freq.o
	$(COMPILE) set_12_gcnmp_vars.f90 

set_initial_profiles.o:	set_initial_profiles.f90 nrtype.o param.o iterdbmd.o \
			 bd_condtn.o numbrs.o ions.o io.o soln.o solcon.o \
	                rfmod.o  consts.o  sourc.o gpsi.o
	$(COMPILE) set_initial_profiles.f90 

set_eqdsk_statefile_values.o:	set_eqdsk_statefile_values.f90 nrtype.o gpsi.o \
				plasma_prop.o cubic_spline_class.o vector_class.o \
				grid_class.o
	$(COMPILE) $< -o $@ 

vector_class.mod:	vector_class.o
vector_class.o:	vector_class.f90 nrtype.o
	$(COMPILE) $< -o $@ 	

tspline.mod:	tspline.o
tspline.o:	tspline.f90 nrtype.o error_handler.o
	$(COMPILE) $< -o $@ 	

zen_state.mod:	zen_state.o
zen_state.o:	zen_state.f90 $(INCLUDES_zen_state)
	$(COMPILE) $< -o $@ 	


#escZ.o:	escZ.c
#	gcc -c escZ.c
#
#2esc.o:	2esc.c
#	gcc -c 2esc.c

#----------------------------------------------------------------
#NUBEAM RELATED ITEMS
# note that the *.o versions of these must be in the libraries that
# are searched at link time.
#----------------------------------------------------------

nubeam_driver.o:    nubeam_driver.f90  $(INCLUDES_nubeam_driver)
	$(COMPILE) nubeam_driver.f90  

nubeam_sub.o:    nubeam_sub.f90 
	$(COMPILE) nubeam_sub.f90

dnubeam_mod.o:    dnubeam_mod.F90 nrtype.o io.o
	$(COMPILE) -DONETWO  $< -o $@ 

dnubeam_driver.o:    dnubeam_driver.f90  $(INCLUDES_nubeam_driver)
	$(COMPILE) dnubeam_driver.f90  

periodic_table_mod.o:    periodic_table_mod.f90 
	$(COMPILE_ntcc) periodic_table_mod.f90

nfreya_namelist.o:	nfreya_namelist.F90 $(INCLUDES_P_nfreya)
	$(COMPILE) -DONETWO  $< -o $@

Nfreya_output.o:	Nfreya_output.f90 transp.o nrtype.o nub2.o nub.o \
			P_Nfreya_12_interface.o  io.o mesh.o geom.o consts.o \
		        plasma_prop.o numbrs.o param.o
	$(COMPILE) Nfreya_output.f90


P_Nfreya_12_interface.o:	P_Nfreya_12_interface.f90 nub.o nub2.o nrtype.o \
			        transp.o io_gcnmp.o solcon.o param.o nf_param.o \
				ions.o fusion.o nub3.o P_Nfreya_nubeam_namelist.o \
                                neutral_beams.o
	$(COMPILE) -DONETWO  $< -o $@

P_Nfreya_nubeam_namelist.o:	P_Nfreya_nubeam_namelist.F90  nrtype.o \
			        transp.o io_gcnmp.o solcon.o param.o nf_param.o \
				ions.o fusion.o nub3.o	dnubeam_mod.o
	$(COMPILE)  $< -o $@

mmm_in.o:	mmm_in.F90 nrtype.o modmmm7_1.o io_gcnmp.o grid_class.o \
		solcon_gcnmp.o error_handler.o consts.o dep_var_gcnmp.o       \
		ions_gcnmp.o neutral_beams.o curden_terms.o plasma_contours.o
	$(COMPILE) -DONETWO   $^

curden_terms.mod:         curden_terms.o  # used only by mmm_in.o, not used in Onetwo
curden_terms.o:         curden_terms.F90   nrtype.o  consts.o grid_class.o
	$(COMPILE) $< -o $@ 

modmmm7_1.o:	 modmmm7_1.f90 nrtype.o w20mod.o
	$(COMPILE)  $< -o $@

w20mod.o:	w20mod.f90 nrtype.o
	$(COMPILE)  $< -o $@

nf_param.o:	nf_param.f90
	$(COMPILE) $< -o $@

nubeam_12.o:	nubeam_12.f $(INCLUDES_nubeam_12)
ifeq ($(NUBEAM_12_BUILD),FALSE)
	$(COMPILE) $(NUBEAM12)
	mv $(basename $(NUBEAM12)).o  nubeam_12.o
else 
	$(COMPILE) nubeam_12.f 
endif

line_av_den.o:	line_av_den.F90 nrtype.o bicubic_spline.o plasma_prop.o \
		consts.o grid_class.o cubic_spline_class.o vector_class.o
	$(COMPILE) $(GFLAGS) $<

#forcebal gcnmp interface routines:
set_forcebal_input.o:	set_forcebal_input.F90 nrtype.o solcon_gcnmp.o \
			consts.o neutral_beams.o plasma_prop.o         \
			grid_class.o MPI_data.o solcon_gcnmp.o
	$(COMPILE) -DONETWO   $<
mmm_input_defaults.o:	mmm_input_defaults.F90 nrtype.o mmm_in.o consts.o \
			modmmm7_1.o
	$(COMPILE) $<

load_input.o:         load_input.F90 iterdbmd_gcnmp.o io_gcnmp.o ions_gcnmp.o \
		      solcon_gcnmp.o grid_class.o bc_values_gcnmp.o nrtype.o \
	              error_handler.o MPI_data.o dep_var_gcnmp.o              \
		      gcnmp_namelist.o neutral_beams.o scale_profiles.o consts.o \
		      plasma_shape.o
	$(COMPILE) -DONETWO $<

get_tot_pressure.o:	get_tot_pressure.f90 nrtype.o grid_class.o ions_gcnmp.o fast_ion_data_gcnmp.o \
			consts.o dep_var_gcnmp.o neutral_beams.o error_handler.o io_gcnmp.o \
			vector_class.o
	$(COMPILE) $< -o $@

scale_profiles.o:	scale_profiles.f90 nrtype.o plasma_prop.o \
			ions_gcnmp.o consts.o grid_class.o  error_handler.o \
		        io_gcnmp.o MPI_data.o 
	$(COMPILE) $<

plasma_shape.o:	plasma_shape.F90 nrtype.o plasma_prop.o consts.o error_handler.o \
		io_gcnmp.o MPI_data.o
	$(COMPILE)  $(GFLAGS) $<

plasma_contours.o:	plasma_contours.f90 nrtype.o plasma_prop.o consts.o error_handler.o \
			io_gcnmp.o bicubic_spline.o grid_class.o
		$(COMPILE)   $< 

#setup_spawn_P_Nfreya.o:		setup_spawn_P_Nfreya.f90 nrtype.o  nfreya_namelist.o
#	$(COMPILE) $< -o $@

P_Nfreya_rpc_interface.o:	P_Nfreya_rpc_interface.f90 nrtype.o  nfreya_namelist.o \
				error_handler.o io_gcnmp.o string_util.o  \
				P_Nfreya_12_interface.o  \
			        ext_prog_info.o iterdbmd_gcnmp.o gcnmp_input.o iterdbmd.o
	$(COMPILE) $< -o $@


ufiles_12.o:    ufiles_12.f90 $(INCLUDES_ufiles_12)
ifeq ($(UBUILD),FALSE)
	$(COMPILE) $(UFILES12)
	@echo ufiles12 a  = $(UFILES12)
	mv $(basename $(UFILES12)).o  ufiles_12.o
else
ifeq ($(FORT),lf95)
	$(COMPILELF95) ufiles_12.f90 
else
	$(COMPILE) ufiles_12.f90 
	@echo ufiles12 b = $(UFILES12)
endif
endif

nbnamelist.mod:	 nbnamelist.o 
nbnamelist.o:    nbnamelist.F90  nrtype.o param.o transp.o numbrs.o
	$(COMPILE)  -DONETWO  $< -o $@  

nbi_restart.mod:	 nbi_restart.o
nbi_restart.o:    nbi_restart.f90 $(INCLUDES_nbi_restart)
ifeq ($(MACHINE),taurus)
	$(COMPILENODB) nbi_restart.f90
else
	$(COMPILE) nbi_restart.f90  
	#mv $(basename $(nbi_restart)).o  nbi_restart.o
endif

get_nboutput.o:    get_nboutput.f90 $(INCLUDES_get_nboutput)
	$(COMPILE) get_nboutput.f90 

transp.mod:  transp.o 
transp.o:  transp.f90  beam_structure.o fast_ion_diffusion.o
	$(COMPILE_nubeam) transp.f90


events.mod:  events.o 
events.o:  events.f90  param.o file_proc.o rfmod.o mhdcom.o tdem.o mhdpar.o
	$(COMPILE_nubeam) events.f90
ifeq ($(HOST),HP-UX)
	cp EVENTS.mod events.mod
endif

beam_structure.mod:  beam_structure.o 
beam_structure.o:  beam_structure.F90  allocate_err.o
	#$(COMPILE_nubeam) beam_structure.F90
	$(COMPILE_nubeam) -DONETWO    $< -o  $@


ifeq ($(HOST),HP-UX)
	cp BEAM_STRUCTURE.mod beam_structure.mod
endif


#modules from transp(NOT from nubeam)
xpload_module.mod:  xpload_module.o
xpload_module.o:  xpload_module.f90  
	$(COMPILE_nubeam) xpload_module.f90

#end nubeam modules --------------------------------------------------------




$(LIBS)/method_lines_library.a: $(LIBS)/method_lines_library.o
	$(DELETE)  $(LIBS)/method_lines_library.a
	$(ARCHIVE) $(LIBS)/method_lines_library.a   $(LIBS)/method_lines_library.o
$(LIBS)/method_lines_library.o: $(LIBS)/method_lines_library.f
	cd         $(LIBS) ; \
	$(COMPILE) $(LIBS)/method_lines_library.f




$(LIBS)/cortlib.a: $(LIBS)/cortlib.o
	$(DELETE)  $(LIBS)/cortlib.a
	$(ARCHIVE) $(LIBS)/cortlib.a   $(LIBS)/cortlib.o
$(LIBS)/cortlib.o: $(LIBS)/cortlib.c
	cd         $(LIBS) ; \
	$(PGCC)  -c    $(LIBS)/cortlib.c




$(LIBS)/$(PORTLIB).a: $(LIBS)/$(PORTLIB).o
	$(DELETE)  $(LIBS)/$(PORTLIB).a
	$(ARCHIVE) $(LIBS)/$(PORTLIB).a   $(LIBS)/$(PORTLIB).o
$(LIBS)/$(PORTLIB).o: $(LIBS)/$(PORTLIB).F
	cd         $(LIBS) ; \
	$(COMPILE) -I../ $(LIBS)/$(PORTLIB).F



##$(LIBS)/$(IMSL).a:    $(LIBS)/$(IMSL).o
##	$(DELETE)  $(LIBS)/$(IMSL).a
##	$(ARCHIVE) $(LIBS)/$(IMSL).a      $(LIBS)/$(IMSL).o
##$(LIBS)/$(IMSL).o:    $(LIBS)/$(IMSL).f
##	cd         $(LIBS) ; \
##	$(COMPILE) $(LIBS)/$(IMSL).f


.PHONY:	clean
clean:
	@echo cleaning up files on $(MACHINE)
	 $(DELETE) $(NAME) *.a *.o *.mod
	 $(DELETE) $(LIBS)/*.o $(LIBS)/*.a 
	 $(DELETE) $(SHARED_MODULES)/*.o
	 #@echo $(LIBS)

data_files: 
	cd $(ONETWO)
