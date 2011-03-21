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


libdir = $(prefix)/lib
noinst_LTLIBRARIES = libgacodemath.la
libgacodemath_la_SOURCES = \
  math_constants.f90 \
  j0y0.f \
  i0.f \
  erf.f \
  p32.f90 \
  invert_p32.f90 \
  gauss_legendre.f90 \
  gauss_integ.f90 \
  polydiff.f90 \
  poly2diff.f90 \
  pascal.f90 \
  bound_deriv.f90 \
  bound_extrap.f90 \
  remap_grid.f90 \
  rjbesl.f \
  ribesl.f \
  cub_spline.f90 \
  cub_spline_deriv.f90 \
  cub_spline_deriv2.f90 \
  gamma.f \
  energy_integral.f90 \
  curtis_rad.f90

if HAVE_MPI
libgacodemath_la_SOURCES += proc_time.f90
endif 
	

# To get to work with optimization reducing compiler wrappers
# or simply for directories with mixed F77 and FC
AM_LIBTOOLFLAGS = --tag=FC
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG)
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG)
CLEANFILES = *.mod *.o *.lo *.lst *.la
