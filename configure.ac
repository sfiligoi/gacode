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

AC_INIT(gacodes, 8.3, srinath@txcorp.com)
AC_PREREQ([2.65])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([config])
dnl AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE
AC_CONFIG_MACRO_DIR([config])

AC_PREFIX_DEFAULT(/internal/gacodes) dnl Install in "/usr/local"

dnl ######################################################################
dnl
dnl Basic initializations
dnl
dnl ######################################################################

builtin(include, config/tx_init.m4)
#AC_DEFINE_UNQUOTED([GACODES_URL], "$PROJECT_URL", "GA_TRANSPORT Project URL")
#AC_DEFINE_UNQUOTED([GACODES_REV], "$PROJECT_REV", "GA_TRANSPORT Revision")

dnl ######################################################################
dnl
dnl Compilers and flags
dnl
dnl ######################################################################

NO_HDF5=true; parallel="yes"
DOUBLE="yes"
builtin(include, config/parcompslibs.m4)
# jrc: why the below is specific to fmcfm?
FCFLAGS="$FCFLAGS $FC_DBL_FLAG"
dnl if test x"$parallel" != "xyes"; then
dnl AC_MSG_ERROR([You cannot build the GA transport codes with out MPI.  First try "--enable-parallel".])
dnl fi

######################################################################
#
# Find fortran information
#
######################################################################

builtin(include, config/tx_fc_aux.m4)

dnl-----------------------------------------------------------------------
dnl Checks for typedefs, structures, and compiler characteristics.
dnl-----------------------------------------------------------------------

dnl jrc: is any of this needed?
if false; then
AC_C_CONST
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_C_VOLATILE
AC_HEADER_STDC
AC_HEADER_SYS_WAIT
AC_CHECK_HEADERS([stdlib.h string.h])
AM_CONFIG_HEADER(config.h)
# Checks for library functions.
AC_FUNC_ALLOCA
CFLAGS="$CFLAGS $C_PIC_FLAG"
CXXFLAGS="$CXXFLAGS $C_PIC_FLAG"
fi

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
builtin(include, config/contrib/acx_blas.m4)
ACX_BLAS(AC_MSG_NOTICE(BLAS was found),AC_MSG_ERROR(BLAS not found.  Use --with-blas=<lib> to set library or --with-blas-dir=<dir> to set search directory.))

dnl Look for LAPACK
builtin(include, config/contrib/acx_lapack.m4)
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

builtin(include,config/python.m4)
AM_CONDITIONAL(CONDPYTHON, test x$python = xtrue)

dnl-----------------------------------------------------------------------
dnl Find the I/O libraries
dnl GYRO can use netcdf via gkv
dnl C. Holland uses the iterdb 2 gyro functionality, and it needs NETCDF
dnl NETCDF_PATH=$HOME/software/netcdf:/contrib/netcdf:/usr/local/netcdf:/opt/netcdf
dnl NETCDF_SP=$HOME:$HOME/$UNIXFLAVOR:$HOME/software:$SUPRA_SEARCH_PATH
dnl-----------------------------------------------------------------------

USE_IO_INFO="true"
builtin(include, config/tx_fciowrappers.m4)
if test "x$NETCDF_LIBS" = "x"; then
   AC_MSG_ERROR([netcdf was not found. NOT SURE HOW TO SET NETCDF LOCATION! (TO FIX)])
fi
# CPPF90_FLAGS="$FPPFLAGS $CPPFLAGS"
CPPF90_FLAGS="$CPP_HAVE_MPI $FLUSH_MACRO $FPPFLAGS"
AC_SUBST(CPPF90_FLAGS)
dnl ######################################################################
dnl
dnl Find documentation executables
dnl
dnl ######################################################################

builtin(include, config/tx_doxygen.m4)
builtin(include, config/tx_latex.m4)

dnl-----------------------------------------------------------------------
dnl Each subdirectory has its own Makefile and dependency
dnl-----------------------------------------------------------------------
AC_CONFIG_FILES([Makefile
    config/Makefile
    config/contrib/Makefile
    envvars/Makefile
    envvars/ga_transport.sh
    envvars/ga_transport.csh
    shared/Makefile
    shared/math/Makefile
    shared/GEO/Makefile
    shared/UMFPACK/Makefile
    shared/EXPRO/Makefile
    shared/profiles_gen/Makefile
    shared/profiles_gen/astra/Makefile
    shared/profiles_gen/modinput/Makefile
    shared/profiles_gen/src/Makefile
    shared/profiles_gen/transp/Makefile
    shared/python/Makefile
    shared/bin/Makefile
    shared/exec/Makefile
    shared/install/Makefile
    gato/Makefile
    gato/data/Makefile
    gato/install/Makefile
    gyro/Makefile
    gyro/BLEND/Makefile
    gyro/TRANSPOSE/Makefile
    gyro/SSUB/Makefile
    gyro/src/Makefile
    gyro/install/Makefile
    gyro/peek/Makefile
    gyro/tools/Makefile
    gyro/tools/PostGYRO/Makefile
    gyro/tools/fieldline/Makefile
    gyro/tools/gyro2gkv/Makefile
    gyro/tools/spectra/Makefile
    gyro/tools/spectra/src/Makefile
    gyro/tools/summary/Makefile
    gyro/bin/Makefile
    gyro/vugyro/Makefile
    neo/Makefile
    neo/bin/Makefile
    neo/nclass/Makefile
    neo/src/Makefile
    tglf/Makefile
    tglf/bin/Makefile
])

#    neo/python/Makefile
#    neo/sim/Makefile
#    neo/src/Makefile
#    gyro/tools/input/Makefile
#    gyro/tools/input/eig01/Makefile
#    gyro/tools/input/eig02/Makefile
#    gyro/tools/input/nl01/Makefile
#    gyro/tools/input/reg01/Makefile
#    gyro/tools/input/reg02/Makefile
#    gyro/tools/input/reg03/Makefile
#    gyro/tools/input/reg04/Makefile
#    gyro/tools/input/reg05/Makefile
#    gyro/tools/input/reg06/Makefile
#    gyro/tools/input/reg07/Makefile
#    gyro/tools/input/reg08/Makefile
#    gyro/tools/input/reg09/Makefile
#    gyro/tools/input/reg10/Makefile
#    gyro/tools/input/reg11/Makefile
#    gyro/tools/input/reg12/Makefile
#    gyro/tools/input/reg13/Makefile
#    gyro/tools/input/reg14/Makefile
#    gyro/tools/input/reg15/Makefile
#    gyro/tools/input/reg16/Makefile
#    gyro/tools/input/reg17/Makefile
#    gyro/tools/input/reg18/Makefile
#    gyro/tools/input/reg19/Makefile
#    gyro/tools/iterdb_read/Makefile
#    gyro/tools/iterdb_read/astra/Makefile
#    gyro/tools/iterdb_read/transp/Makefile
#    gyro/tools/iterdb_read/gato/Makefile
#    gyro/tools/iterdb_read/src/Makefile
#    gyro/tools/vtkgyro/Makefile
#    gyro/tools/vtkgyro/common/Makefile
#    gyro/tools/vtkgyro/movie/Makefile
#    gyro/tools/vtkgyro/square/Makefile
#    gyro/tools/vtkgyro/square_2/Makefile
#    gyro/tools/vtkgyro/torus/Makefile
#    gyro/tools/vtkgyro/torus-old/Makefile
#    gyro/tools/vtkgyro/torus-old/sample-INPUT/Makefile
#    gyro/tools/vtkgyro/volume/Makefile
#    gyro/batch/Makefile
#    gyro/exec/Makefile
#    gyro/python/Makefile
#    gyro/install/Makefile
#    gyro/sim/Makefile
#    gyro/sim/test/Makefile
#    gyro/sim/test/reg01/Makefile
#    gyro/sim/test/reg02/Makefile
#    gyro/sim/test/reg03/Makefile
#    gyro/sim/test/reg04/Makefile
#    gyro/sim/test/reg05/Makefile
#    gyro/sim/test/reg06/Makefile
#    gyro/sim/test/reg07/Makefile
#    gyro/sim/test/reg08/Makefile
#    gyro/sim/test/reg09/Makefile
#    gyro/sim/test/reg10/Makefile
#    gyro/sim/test/reg11/Makefile
#    gyro/sim/test/reg12/Makefile
#    gyro/sim/test/reg13/Makefile
#    gyro/sim/test/reg14/Makefile
#    gyro/sim/test/reg15/Makefile
#    gyro/sim/test/reg16/Makefile
#    gyro/sim/test/reg17/Makefile
#    neo/tools/Makefile
#    neo/tools/vgen/Makefile
#    neo/tools/vgen/src/Makefile
#    neo/tools/squeeze/Makefile
#    neo/tools/squeeze/src/Makefile
#    neo/tools/input/Makefile
#    neo/tools/input/ade_c1/Makefile
#    neo/tools/input/ade_c2/Makefile
#    neo/tools/input/ade_c3/Makefile
#    neo/tools/input/ade_miller/Makefile
#    neo/tools/input/exp_geo_norot/Makefile
#    neo/tools/input/exp_geo_rot/Makefile
#    neo/tools/input/exp_geo_rot_multir/Makefile
#    neo/tools/input/exp_miller_norot/Makefile
#    neo/tools/input/exp_miller_rot/Makefile
#    neo/tools/input/fow_ade/Makefile
#    neo/tools/input/fow_exp_ade/Makefile
#    neo/tools/input/fow_exp_imp/Makefile
#    neo/tools/input/gke_c1/Makefile
#    neo/tools/input/gke_c2/Makefile
#    neo/tools/input/gke_c3/Makefile
#    neo/tools/input/imp_c1/Makefile
#    neo/tools/input/imp_c2/Makefile
#    neo/tools/input/imp_c3/Makefile
#    neo/tools/input/rot_ade/Makefile
#    neo/tools/input/rot_ade_shear/Makefile
#    neo/tools/input/vgen/Makefile
#    neo/tools/input/vgen1/Makefile
#    neo/tools/input/vgen_default/Makefile
#    tglf/Makefile
#    tglf/src/Makefile
#    tglf/bin/Makefile
#    tglf/python/Makefile
#    tglf/sim/Makefile
#    tglf/sim/ttest/Makefile
#    tgyro/Makefile
#    tgyro/bin/Makefile
#    tgyro/math/Makefile
#    tgyro/python/Makefile
#    tgyro/qfm/Makefile
#    tgyro/qfm/src/Makefile
#    tgyro/sim/Makefile
#    tgyro/src/Makefile
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

builtin(include, config/rminterlibdeps.ac)

dnl-----------------------------------------------------------------------
dnl  Print out config summary file.
dnl-----------------------------------------------------------------------

TX_AC_SUMMARY_FINALIZE

dnl-----------------------------------------------------------------------
dnl  Output of files must be last
dnl-----------------------------------------------------------------------

AC_OUTPUT

