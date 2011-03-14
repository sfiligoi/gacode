## ##########################################################################
##
## File:		Makefile.am
##
## Author: Mark Fahey
##
## Purpose:  Process this file with automake to produce Makefile
##
## Version:	$Id: Makefile.am 2 2009-09-01 23:34:41Z srinath $
##
## ############################################################################


libdir = $(prefix)/lib
noinst_LTLIBRARIES = libgacodeumfpack.la
HARWELL	= mc21b.f mc13e.f
UMFD = \
	umd21i.f \
	umd2fa.f \
	umd2p1.f \
	umd2rg.f \
	umz21i.f \
	umz2fa.f \
	umz2p2.f \
	umz2s2.f \
	umd2co.f \
	umd2fb.f \
	umd2p2.f \
	umd2s2.f \
	umz2co.f \
	umz2fb.f \
	umz2r0.f \
	umz2sl.f \
	umd2er.f \
	umd2fg.f \
	umd2r0.f \
	umd2sl.f \
	umz2er.f \
	umz2fg.f \
	umz2r2.f \
	umz2so.f \
	umd2f0.f \
	umd2in.f \
	umd2r2.f \
	umd2so.f \
	umz2f0.f \
	umz2in.f \
	umz2ra.f \
	umz2su.f \
	umd2f1.f \
	umd2lt.f \
	umd2ra.f \
	umd2su.f \
	umz2f1.f \
	umz2of.f \
	umz2rf.f \
	umd2f2.f \
	umd2of.f \
	umd2rf.f \
	umd2ut.f \
	umz2f2.f \
	umz2p1.f \
	umz2rg.f
       
libgacodeumfpack_la_SOURCES = $(HARWELL) $(UMFD)

# To get to work with optimization reducing compiler wrappers
# or simply to get .f to work with new automake
AM_LIBTOOLFLAGS = --tag=FC
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG)
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG)
CLEANFILES = *.mod *.o *.lo *.lst *.la
