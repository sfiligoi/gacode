## ##########################################################################
##
## File:		Makefile.am
##
## Author: Mark Fahey
##
## Purpose:  Process this file with automake to produce Makefile
##
## Version:	$Id: Makefile.am 176 2010-07-27 15:25:25Z srinath $
##
## ############################################################################


bin_PROGRAMS     = peek

gyro_includes= -I $(top_builddir)/gyro/src

# To get to work with optimization reducing compiler wrappers
IOINC=@HDF5_INC@ @FCIOWRAPPERS_INC@
IOLIBS=@FCIOWRAPPERS_LIBS@ @HDF5_LIBS@ @NETCDF_LIBS@
AM_LIBTOOLFLAGS = --tag=FC
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(IOINC) $(gyro_includes) $(CPPF90_FLAGS) -I.
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(IOINC) $(gyro_includes) $(CPPF90_FLAGS) -I.

peek_SOURCES = \
	       peek_globals.f90 \
	       read_profile_vugyro.f90 \
	       read_time.f90 \
         peek.f90 

peek_LDADD  = $(top_builddir)/gyro/src/libgacodegyro_SRC.la
#	@FFTW_LIBS@  @LAPACK_LIBS@  @BLAS_LIBS@ $(IOLIBS) $(Z_LIBS) \
#	$(top_builddir)/gyro/EXPRO/libEXPRO.la
#
peek_INCLUDES = -I $(top_builddir)/gyro/src 

if TX_FORTRAN_MODNAMEISCAP
nodist_include_HEADERS = PEEK_GLOBALS.$(TX_FORTRAN_MODEXT)  
PEEK_GLOBALS.$(TX_FORTRAN_MODEXT):   peek_globals.f90
else
nodist_include_HEADERS = peek_globals.$(TX_FORTRAN_MODEXT) 
peek_globals.$(TX_FORTRAN_MODEXT):   peek_globals.f90
endif


CLEANFILES = *.mod *.f90 *.lo *.lst *.la

