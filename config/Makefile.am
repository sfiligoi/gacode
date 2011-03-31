## ##########################################################################
##
## File:	Makefile.am
##
## Purpose:	Run automake at the top level to make Makefile.in in all
##		subdirectories.
##
## Version:	$Id: Makefile.am 3778 2011-01-13 04:33:29Z dws $
##
## Copyright &copy; 1999-2002, University of Colorado
## Copyright &copy; 2002-2010, University of Colorado and Tech-X Corporation
## See LICENSE file for conditions of use.
##
## ##########################################################################

# Following is config.guess config.sub depcomp install-sh missing mkinstalldirs
# along with Makefile.in *.am *.sh *.m4 *.pl

SUBDIRS		= contrib

CONFIG_DIST	= config.guess config.sub depcomp install-sh missing \
			mkinstalldirs

DIST_FILES = COPYING
DIST_SUFFIXES = .am .m4 .pl .sh
NODIST_FILES = Makefile.am
EXTRA_DIST = \
   COPYING \
   depend.ac rminterlibdeps.ac \
   bindist.am bundle.am clean-ac.am clean-cxx.am depend.am \
   doxygen.am filelists.am \
   docs.am sphinx.am \
   ant.m4 aztec.m4 blas.m4 boost.m4 ccafe.m4 ccfits.m4 cfitsio.m4 \
   cmee.m4 corba.m4 cxx.m4 defaultcomps.m4 f77.m4 f77flags.m4 f77interop.m4 \
   fc.m4 fcflags.m4 fftw.m4 finterop.m4 flags.m4 \
   globus.m4 gsoap.m4 hdf5_cpp.m4 hl5.m4 idldepend.m4 ionpack.m4 ittidl.m4 \
   java.m4 java2.m4 jchart.m4 jcorba.m4 jmflibs.m4 jun4java.m4 \
   lapack.m4 libs.m4 linlibs.m4 macros.m4 mdsplus.m4 \
   mpif77.m4 mpilauncher.m4 netcdf_fortran.m4 nm.m4 nsecc.m4 ntcbase.m4 \
   numpy.m4 opendx.m4 opengl.m4 os_tag.m4 parcompslibs.m4 pvwave.m4 pycxx.m4 \
   python.m4 qscimpl.m4 qt3.m4 qt4.m4 rsiidl.m4 serialcomps.m4 sshtools.m4 \
   swig.m4 tx_tau.m4 tcltk.m4 trilinos.m4 tx_ac_summary.m4 tx_babel.m4 \
   tx_blaslapack.m4 tx_check_args.m4 tx_doxygen.m4 \
   tx_fc_allocatable_component.m4 tx_fc_aux.m4 \
   tx_fc_bindc.m4 tx_fc_find_main.m4 tx_fc_has_intrinsic.m4 \
   tx_fc_iso_c_binding.m4 tx_fc_module_capitalization.m4 tx_fc_netcdf.m4 \
   tx_fc_sizeof_array_descriptor.m4 tx_fc_sizeof.m4 \
   tx_fc_sizeof_pointer_container.m4 tx_gpulib.m4 \
   tx_init.m4 tx_latex.m4 tx_llnl_iso_c_binding.m4 tx_netcdf.m4 tx_superlu.m4 \
   tx_svn_info.m4 tx_sz.m4 tx_petsc.m4 tx_hdf5.m4 \
   tx_licmgr.m4 tx_netlib_lite.m4 tx_ntcctransport.m4 tx_fciowrappers.m4 tx_transplibs.m4 \
   tx_uedge.m4 tx_visit.m4 txbase.m4 txphysics.m4 txqplot.m4 txsearch.m4 \
   vorpal.m4 \
   waspc.m4 wx.m4 x11.m4 x11_ext.m4 x11xm.m4 xerces.m4 xpm.m4 \
   cvs2cl.pl fdepend.pl writedep.pl \
   rmcomms.sed \
   automake.sh cleanconf.sh \
   mpifuncs.sh regenconf.sh updatemfs.sh

noinst_HEADERS =

cleanautoconf:
	rm -f Makefile.in
	if test -n "$(SUBDIRS)"; then \
	  for i in $(SUBDIRS); do cd $$i; $(MAKE) cleanautoconf; done \
	fi
	rm -f Makefile.bak
	mv Makefile Makefile.bak

.PHONY: cleanautoconf

## #######
##
## Other targets
##
## #######

include filelists.am
