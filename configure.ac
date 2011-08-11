dnl                                               -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.
dnl -----------------------------------------------------------------------
dnl
dnl File:    configure.in
dnl
dnl Purpose:  configuration script for GA's gyro
dnl
dnl Copyright 2005-2010, Tech-X Corporation, All Rights Reserved.
dnl
dnl Version:  $Id: configure.ac 1267 2009-07-28 13:10:59Z srinath$
dnl
dnl------------------------------------------------------------------------

AC_INIT(gacodes, 8.3.8, srinath@txcorp.com)
AC_PREREQ([2.65])
AC_CONFIG_AUX_DIR(shared/autotool_config)
AC_CONFIG_MACRO_DIR([shared/autotool_config])
dnl AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE
AC_CONFIG_MACRO_DIR([shared/autotool_config])

AC_PREFIX_DEFAULT(/internal/gacodes) dnl Install in "/usr/local"

dnl ######################################################################
dnl
dnl Basic initializations
dnl
dnl ######################################################################

builtin(include, shared/autotool_config/tx_init.m4)
#AC_DEFINE_UNQUOTED([GACODES_URL], "$PROJECT_URL", "GA_TRANSPORT Project URL")
#AC_DEFINE_UNQUOTED([GACODES_REV], "$PROJECT_REV", "GA_TRANSPORT Revision")

dnl ######################################################################
dnl
dnl Compilers and flags
dnl
dnl ######################################################################

NO_HDF5=true; parallel="yes"
DOUBLE="yes"
builtin(include, shared/autotool_config/parcompslibs.m4)


if test x"$parallel" = "xyes"; then
  FCPREFLAG='-D__MPI'
fi
FCFLAGS="$FCFLAGS $FC_DBL_FLAG $FCPREFLAG"

######################################################################
#
# Find fortran information
#
######################################################################

builtin(include, shared/autotool_config/tx_fc_aux.m4)

dnl ######################################################################
dnl
dnl Libtool (may get additional flags)
dnl
dnl ######################################################################

LT_INIT
AC_PROG_LIBTOOL

dnl ######################################################################
dnl
dnl Extra flags
dnl
dnl ######################################################################



dnl J.Candy recommends not to build any code in Eigen (with PETSC dependence)
dnl and any code in the tools directory.  This means everything will be moved to EXTRA_DIST

dnl-----------------------------------------------------------------------
dnl Math libraries
dnl-----------------------------------------------------------------------
echo  >> $config_summary_file

echo  "Linear algebra libraries used by GYRO:" >> $config_summary_file
builtin(include, shared/autotool_config/contrib/acx_blas.m4)
ACX_BLAS(AC_MSG_NOTICE(BLAS was found),AC_MSG_ERROR(BLAS not found.  Use --with-blas=<lib> to set library or --with-blas-dir=<dir> to set search directory.))

dnl Look for LAPACK
builtin(include, shared/autotool_config/contrib/acx_lapack.m4)
ACX_LAPACK(AC_MSG_NOTICE(LAPACK was found),
AC_MSG_ERROR(LAPACK not found.  Use --with-lapack=<lib> to set library or --with-lapack-dir=<lapackdir> to set search directory.))

unset MUMPS_PATH
for i in `echo $SUPRA_SEARCH_PATH | tr ':' ' '`; do
  MUMPS_PATH="$MUMPS_PATH:$i/mumps"
done
TX_LOCATE_PKG( [MUMPS], [$MUMPS_PATH], [mumps.h], [mumps])

unset ESSL_PATH
for i in `echo $SUPRA_SEARCH_PATH | tr ':' ' '`; do
  ESSL_PATH="$ESSL_PATH:$i/essl"
done
TX_LOCATE_PKG( [ESSL], [$ESSL_PATH], [essl.h], [essl])

unset LIBSCI_PATH
for i in `echo $SUPRA_SEARCH_PATH | tr ':' ' '`; do
  LIBSCI_PATH="$LIBSCI_PATH:$i/libsci"
done
TX_LOCATE_PKG( [LIBSCI], [$LIBSCI_PATH], [libsci.h], [libsci])

dnl this is added for mpiio possible usauge
AM_CONDITIONAL([HAVE_MPIIO],[false])
AM_CONDITIONAL([HAVE_XE1],[false])

dnl ######################################################################
dnl  fftw is optional for gyro
dnl ######################################################################

unset FFTW_PATH
for i in `echo $SUPRA_SEARCH_PATH | tr ':' ' '`; do
  FFTW_PATH="$FFTW_PATH:$i/fftw"
done
TX_LOCATE_PKG(
  [FFTW],
  [$FFTW_PATH],
  [fftw.h,rfftw.h],
  [rfftw,fftw])


dnl ######################################################################
dnl
dnl Python
dnl
dnl ######################################################################

builtin(include,shared/autotool_config/python.m4)
AM_CONDITIONAL(CONDPYTHON, test x$python = xtrue)

dnl-----------------------------------------------------------------------
dnl Find the I/O libraries
dnl GYRO can use netcdf via gkv
dnl C. Holland uses the iterdb 2 gyro functionality, and it needs NETCDF
dnl NETCDF_PATH=$HOME/software/netcdf:/contrib/netcdf:/usr/local/netcdf:/opt/netcdf
dnl NETCDF_SP=$HOME:$HOME/$UNIXFLAVOR:$HOME/software:$SUPRA_SEARCH_PATH
dnl-----------------------------------------------------------------------

USE_IO_INFO="true"
ac_use_fortran_hdf5="yes"
builtin(include, shared/autotool_config/tx_fciowrappers.m4)
builtin(include, shared/autotool_config/tx_hdf5.m4)
builtin(include, shared/autotool_config/tx_netcdf.m4)

AC_ARG_WITH(fciowrappers,
 AC_HELP_STRING([--with-fciowrappers], [to build with fciowrappers]),
        usefciowrappers=$withval, usefciowrappers=no)
AM_CONDITIONAL(USE_FCIOWRAPPERS, test x$usefciowrappers = xyes )

if test x$usefciowrappers = xyes;then
  AC_DEFINE(USE_FCIOWRAPPERS, [], [Defined if fciowrappers to be used.])
fi


CPPF90_FLAGS="$CPP_HAVE_MPI $FLUSH_MACRO $FPPFLAGS"
AC_SUBST(CPPF90_FLAGS)
dnl ######################################################################
dnl
dnl Find documentation executables
dnl
dnl ######################################################################

builtin(include, shared/autotool_config/tx_doxygen.m4)
builtin(include, shared/autotool_config/tx_latex.m4)

dnl-----------------------------------------------------------------------
dnl Each subdirectory has its own Makefile and dependency
dnl-----------------------------------------------------------------------
AC_CONFIG_FILES([Makefile
    shared/autotool_config/Makefile
    shared/autotool_config/contrib/Makefile
    shared/envvars/Makefile
    shared/envvars/ga_transport.sh
    shared/envvars/ga_transport.csh
    shared/Makefile
    shared/math/Makefile
    shared/GEO/Makefile
    shared/UMFPACK/Makefile
    shared/EXPRO/Makefile
    shared/vshdf5/Makefile
    shared/profiles_gen/Makefile
    shared/profiles_gen/astra/Makefile
    shared/profiles_gen/modinput/Makefile
    shared/profiles_gen/src/Makefile
    shared/profiles_gen/transp/Makefile
    shared/python/Makefile
    shared/bin/Makefile
    shared/exec/Makefile
    gyro/Makefile
    gyro/BLEND/Makefile
    gyro/TRANSPOSE/Makefile
    gyro/SSUB/Makefile
    gyro/GKEIGEN/Makefile
    gyro/vshdf5/Makefile
    gyro/src/Makefile
    gyro/peek/Makefile
    gyro/tools/Makefile
    gyro/tools/PostGYRO/Makefile
    gyro/tools/fieldline/Makefile
    gyro/tools/gyro2gkv/Makefile
    gyro/tools/spectra/Makefile
    gyro/tools/spectra/src/Makefile
    gyro/tools/summary/Makefile
    gyro/tools/PlotGYROpy/Makefile
    gyro/bin/Makefile
    gyro/vugyro/Makefile
    neo/Makefile
    neo/bin/Makefile
    neo/nclass/Makefile
    neo/src/Makefile
    tglf/Makefile
    tglf/bin/Makefile
    tglf/python/Makefile
    tglf/src/Makefile
    tgyro/Makefile
    tgyro/bin/Makefile
    tgyro/math/Makefile
    tgyro/qfm/Makefile
    tgyro/src/Makefile
])

#removed from build system
#    gato/Makefile
#    gato/data/Makefile
#    xptor/Makefile
#    xptor/datread/Makefile
#    xptor/misc/Makefile
#    xptor/inc/Makefile
#    xptor/models/Makefile
#    xptor/nclmisc/Makefile
#    xptor/ptmisc/Makefile
#    xptor/solver/Makefile
#    xptor/source/Makefile
#    xptor/driver/Makefile


#    xptor/idl/Makefile
#
#    tgyro/tools/Makefile
#    tgyro/tools/asymptote/Makefile
#    tgyro/tools/asymptote/src/Makefile
#    tgyro/tools/idl/Makefile
#    tgyro/tools/input/Makefile
#    tgyro/tools/input/stest_1/Makefile
#    tgyro/tools/input/stest_2/Makefile
#    tgyro/tools/input/stest_3/Makefile
#    tgyro/tools/input/stest_4/Makefile
#    tgyro/tools/input/stest_5/Makefile
#    tgyro/tools/input/stest_5/GYRO1/Makefile
#    tgyro/tools/input/stest_5/GYRO2/Makefile
#    tgyro/tools/input/stest_5/GYRO3/Makefile
#    tgyro/tools/input/stest_5/GYRO4/Makefile
#    tgyro/tools/input/wtest_1/Makefile
#    tgyro/tools/input/wtest_2/Makefile
#    tgyro/tools/input/wtest_3/Makefile
#    tgyro/tools/input/astra_dz18_tglf/Makefile
#    tgyro/tools/input/d3_20_128913_tglf/Makefile
#    tgyro/tools/input/d3_4_gyro/Makefile
#    tgyro/tools/input/d3_4_gyro/GYRO1/Makefile
#    tgyro/tools/input/d3_4_gyro/GYRO2/Makefile
#    tgyro/tools/input/d3_4_gyro/GYRO3/Makefile
#    tgyro/tools/input/d3_4_gyro/GYRO4/Makefile
#    tgyro/tools/input/d3_4_tglf/Makefile
#    tgyro/tools/input/d3_8_128913_tglfneo/Makefile
#    ])
dnl ######################################################################
dnl
dnl  Remove interlibrary dependency checking
dnl
dnl ######################################################################

builtin(include, shared/autotool_config/rminterlibdeps.ac)

dnl-----------------------------------------------------------------------
dnl  Print out config summary file.
dnl-----------------------------------------------------------------------

TX_AC_SUMMARY_FINALIZE

dnl-----------------------------------------------------------------------
dnl  Output of files must be last
dnl-----------------------------------------------------------------------

AC_OUTPUT

