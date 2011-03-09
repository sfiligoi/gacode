dnl ######################################################################
dnl
dnl File:	tx_hdf5.m4
dnl
dnl Purpose:	Determine the locations of hdf5 includes and libraries
dnl
dnl Version: $Id: tx_hdf5.m4 3804 2011-03-02 23:46:26Z cary $
dnl
dnl Tech-X configure system
dnl
dnl need to bring this into compliance:
dnl	HDF5_LIBDIR is the directory containing -lhdf5
dnl     HDF5_LIBS should be the entire variable needed, i.e.,
dnl		-L$(HDF5_LIBDIR) -lhdf5
dnl	HDF5_LIB set to -lhdf5 if found.
dnl	same should be for INCDIR and INCS.
dnl
dnl This m4 file defines:
dnl    HAVE_HDF5
dnl    HDF5_INCDIR
dnl    HDF5_MODDIR
dnl    HDF5_INC
dnl    HDF5_LIBS
dnl    HDF5_LIBDIR
dnl
dnl  2005-02-30 - SEK
dnl    Added hdf5 enable option for codes that do not require hdf5
dnl ######################################################################

dnl Includes functions from txsearch.m4
builtin(include, config/txsearch.m4)

dnl ######################################################################
dnl
dnl Allow the user to specify an overall hdf5 directory.  If specified,
dnl we look for include and lib under this.
dnl Notes on Fortran support.
dnl  Problem is that hdf5 f90 wrappers need to be compiled with same
dnl    compiler as code compiler.  Perform check using hdf5 settings file
dnl  To use this functionality, one needs to set in configure.ac:
dnl    ac_check_fortran_hdf5=true        # Perform fortran hdf5 tests
dnl
dnl ######################################################################

AC_ARG_WITH(hdf5-dir, AC_HELP_STRING([--with-hdf5-dir],
        [hdf5 installation directory]),
        HDF5_DIR="$withval")

dnl ######################################################################
dnl
dnl Find hdf5 includes - looking in include location if present,
dnl otherwise in dir/include if present, otherwise in default locations,
dnl first parallel, then serial.
dnl
dnl ######################################################################

AC_ARG_WITH(hdf5-incdir, AC_HELP_STRING([--with-hdf5-incdir],
        [hdf5 include directory]),
        HDF5_INCDIR="$withval")

if test -n "$HDF5_INCDIR"; then
  HDF5_PATH=`(cd $HDF5_INCDIR/..; pwd -P)`
elif test -n "$HDF5_DIR"; then
  HDF5_PATH=$HDF5_DIR
else
  HDF5_SP=$SUPRA_SEARCH_PATH
  unset HDF5_SERPATH
  for i in `echo $HDF5_SP | tr ':' ' '`; do
    HDF5_SERPATH="$HDF5_SERPATH:$i/hdf5"
  done
  for i in `echo $HDF5_SP | tr ':' ' '`; do
    HDF5_SERPATH="$HDF5_SERPATH:$i"
  done
  unset HDF5_PARPATH
  for i in `echo $HDF5_SP | tr ':' ' '`; do
    HDF5_PARPATH="$HDF5_PARPATH:$i/hdf5-par:$i/hdf5mpi" # Last is deprecated
  done
  for i in `echo $HDF5_SP | tr ':' ' '`; do
    HDF5_PARPATH="$HDF5_PARPATH:$i" # Last is deprecated
  done
  # echo PAR_HDF5_IS_DEFAULT = $PAR_HDF5_IS_DEFAULT
  if test -z "$PAR_HDF5_IS_DEFAULT"; then
    if test -n "$parallel"; then
      PAR_HDF5_IS_DEFAULT=$parallel
    elif test -n "$paralleldefault"; then
      PAR_HDF5_IS_DEFAULT=$paralleldefault
    else
      PAR_HDF5_IS_DEFAULT=yes
    fi
  else
    PAR_HDF5_IS_DEFAULT=yes
  fi
  if test "$PAR_HDF5_IS_DEFAULT" = yes; then
    HDF5_PATH=$HDF5_PATH:$HDF5_PARPATH:$HDF5_SERPATH
  else
    HDF5_PATH=$HDF5_PATH:$HDF5_SERPATH:$HDF5_PARPATH
  fi
fi

dnl ######################################################################
dnl
dnl Search for includes and libraries.
dnl
dnl ######################################################################

if test "$ac_use_fortran_hdf5" = yes; then
  echo Searching for hdf5 in $HDF5_PATH.
  TX_LOCATE_PKG(
    [HDF5],
    [$HDF5_PATH],
    [hdf5.h],
    [hdf5_fortran, hdf5_f90cstub, hdf5_hl, hdf5],
    [include],
    [lib])
dnl The following line make HAVE_HDF5 true if an hdf5 is found
  AM_CONDITIONAL(HAVE_HDF5, test -n "$HDF5_LIBS")
  if test -n "$HDF5_LIBS"; then
    AC_DEFINE(HAVE_HDF5, [], "Define if HDF5 is found")
  fi
dnl the following search for the presence of lib hdf5_f90cstub
dnl not finding this library is not an error, thus we use a separate search
dnl JRC: putting above for exactly this reason.  It is not an error.
dnl  TX_LOCATE_PKG(
dnl    [HDF5XTRA],
dnl    [$HDF5_PATH],
dnl    [],
dnl    [hdf5_f90cstub],
dnl    [],
dnl    [lib])
  HDF5_LIB="-lhdf5_fortran"
  if test -n "$HDF5_LIB_HDF5_F90CSTUB"; then
    HDF5_LIB="$HDF5_LIB -lhdf5_f90cstub"
  fi
  HDF5_LIB="$HDF5_LIB -lhdf5"
# Look for hdf5.mod or HDF5.mod
# As of hdf5-1.8.5, this file is under the fortran subdir
  HDF5_FC_INC="$HDF5_INC"
  echo HDF5_FC_INC = $HDF5_FC_INC
  TX_PATH_FILES(HDF5_MOD, hdf5.mod HDF5.mod, "", $HDF5_INCDIR/fortran:$HDF5_INCDIR/../fortran:$HDF5_INCDIR:$HDF5_LIBDIR)
  if test -n "$HDF5_MOD"; then
    HDF5_MODDIR=`dirname $HDF5_MOD`
    AC_SUBST(HDF5_MODDIR)
    if test $HDF5_MODDIR != $HDF5_INCDIR; then
      HDF5_FC_INC="$HDF5_FC_INC -I$HDF5_MODDIR"
    fi
    AC_DEFINE(HAVE_HDF5_MOD, [], "Define if have hdf5.mod file")
  else
    echo "Neither of hdf5.mod HDF5.mod found in $HDF5_INCDIR/fortran:$HDF5_INCDIR/../fortran:$HDF5_INCDIR:$HDF5_LIBDIR"
  fi
  echo HDF5_FC_INC = $HDF5_FC_INC
  TX_PRINT_VAR(HDF5_MOD)
  TX_PRINT_VAR(HDF5_MODDIR)
  TX_PRINT_VAR(HDF5_FC_INC)
else
  TX_LOCATE_PKG(
    [HDF5],
    [$HDF5_PATH],
    [hdf5.h],
    [hdf5_hl, hdf5],
    [include],
    [lib])
  HDF5_LIB="-lhdf5"
fi
AC_SUBST(HDF5_LIB)
TX_PRINT_VAR(HDF5_LIB)
AC_SUBST(HDF5_MOD)
AC_SUBST(HDF5_FC_INC)
AM_CONDITIONAL(HAVE_HDF5_MOD, test -n "$HDF5_MOD")

dnl error conditions
if test -z "$HDF5_INC_HDF5_H"; then
  AC_MSG_WARN(hdf5.h not found in $HDF5_PATH/include.  Use --with-hdf5-dir= to specifiy the installation directory.)
  HDF5_INC=" "
  ac_cv_have_hdf5=no
  HDF5_LIB=""
else
  ac_cv_have_hdf5=yes
fi

dnl ##############################################################
dnl
dnl  Enable fortran consistency checks done later
dnl
dnl ##############################################################

dnl JRC: The code below seems to break with each release of HDF5,
dnl so now have a new variable, ac_check_fortran_hdf5, that invokes it.
dnl Some may not want this.
if test $ac_cv_have_hdf5 = yes; then
if test "$ac_require_fortran_hdf5" = "no"; then
    ac_require_hdf5="no"
fi
if test "$ac_check_fortran_hdf5" = yes; then
  HDF5_LIB="hdf5_fortran hdf5"
  if test -f $HDF5_LIBDIR/libhdf5_fortran.a; then
dnl test to make sure hdf5_fortran exists
    if test -f $HDF5_LIBDIR/libhdf5.settings; then
      hdf5_f90_full=`grep 'Fortran Compiler' $HDF5_LIBDIR/libhdf5.settings | cut -f2 -d:`
    fi
    if test -z "$hdf5_f90_full"; then
      # echo Fortran Compiler not found in $HDF5_LIBDIR/libhdf5.settings
      if test -f $HDF5_LIBDIR/libhdf5_fortran.settings; then
        hdf5_f90_full=`grep 'Fortran Compiler' $HDF5_LIBDIR/libhdf5_fortran.settings | sed 's/^.*: *//' | sed -n '1p'`
        # echo hdf5_f90_full = $hdf5_f90_full
      else
        if test "$ac_require_fortran_hdf5" = "no"; then
	   HDF5_LIBDIR=""
	   HDF5_LIBS=""
	   HDF5_INCDIR=""
	   HDF5_MODDIR=""
	   HDF5_INC=""
	   HDF5_FC=""
	else
          AC_MSG_ERROR(Cannot determine Fortran compiler for hdf5.)
	fi
      fi
    fi
# Remove leading spaces and then all after first space
    hdf5_f90_full=`echo $hdf5_f90_full | sed -e 's/^  *//' -e 's/ .*$//'`
    # echo hdf5_f90_full = $hdf5_f90_full
    if test "$ac_require_fortran_hdf5" = "no"; then
      :
    else
     if test -z "$hdf5_f90_full"; then
      AC_MSG_ERROR(Failed to get hdf5_f90_full from hdf5 settings files.)
     fi
    fi
    HDF5_FC=`basename $hdf5_f90_full`
dnl jrc: This is not working, so will quit if bag
    if test -z "$HDF5_FC"; then
        if test "$ac_require_fortran_hdf5" = "no"; then
	   HDF5_LIBDIR=""
	   HDF5_LIBS=""
	   HDF5_INCDIR=""
	   HDF5_MODDIR=""
	   HDF5_INC=""
	   HDF5_FC=""
	  else
         AC_MSG_ERROR(Failed to get HDF5_FC.)
	  fi
    fi
  fi
fi
fi
AC_SUBST(HDF5_FC)

dnl ######################################################################
dnl
dnl See if built parallel
dnl
dnl ######################################################################

if test $ac_cv_have_hdf5 = yes; then
  if test -f $HDF5_INCDIR/H5config.h; then
    hdf5par=`grep "HAVE_PARALLEL 1" $HDF5_INCDIR/H5config.h`
  elif test -f $HDF5_INCDIR/H5pubconf.h; then
    hdf5par=`grep "HAVE_PARALLEL 1" $HDF5_INCDIR/H5pubconf.h`
  fi
fi

dnl ######################################################################
dnl
dnl If built parallel, may need to include parallel file system library
dnl
dnl ######################################################################

if test $ac_cv_have_hdf5 = yes; then
if test -n "$hdf5par"; then
  GPFS_LIBPATH=/usr/lib:/usr/local/lib:/usr/lpp/mmfs/lib
  TX_PATH_FILE(GPFS_LIB, libgpfs.a, "", $GPFS_LIBPATH)
  if test -z "$GPFS_LIB"; then
    TX_PATH_FILE(GPFS_LIB, libgpfs.so, "", $GPFS_LIBPATH)
  fi
# JRC 11Jul08: Should not have?
  if test -n "$GPFS_LIB"; then
    GPFS_LIB=""
    GPFS_LIBDIR=""
    GPFS_LIBS=""
  fi
fi
fi
AC_SUBST(GPFS_LIBS)
AC_SUBST(GPFS_LIB)
AC_SUBST(GPFS_LIBDIR)
TX_PRINT_VAR(GPFS_LIBS)
TX_PRINT_VAR(GPFS_LIB)
TX_PRINT_VAR(GPFS_LIBDIR)

dnl ######################################################################
dnl
dnl Determine version and corresponding changes to interfaces.
dnl
dnl ######################################################################

if test $ac_cv_have_hdf5 = yes; then
  H5_VERS_MAJOR=`grep '#define H5_VERS_MAJOR' $HDF5_INCDIR/H5public.h | sed 's/^.*H5_VERS_MAJOR//' | sed 's/\/.*$//' | tr -d [:space:]`
  H5_VERS_MINOR=`grep '#define H5_VERS_MINOR' $HDF5_INCDIR/H5public.h | sed 's/^.*H5_VERS_MINOR//' | sed 's/\/.*$//' | tr -d [:space:]`
  H5_VERS_RELEASE=`grep '#define H5_VERS_RELEASE' $HDF5_INCDIR/H5public.h | sed 's/^.*H5_VERS_RELEASE//' | sed 's/\/.*$//' | tr -d [:space:]`
  unset new_H5Sselect_hyperslab_ifc
  unset h5_use_16
  if test $H5_VERS_MAJOR -gt 1; then
    new_H5Sselect_hyperslab_ifc=true
  elif test $H5_VERS_MAJOR = 1; then
    if test $H5_VERS_MINOR -gt 6; then
      new_H5Sselect_hyperslab_ifc=true
    elif test $H5_VERS_MINOR = 6; then
      if test $H5_VERS_RELEASE -ge 4; then
        new_H5Sselect_hyperslab_ifc=true
      fi
    fi
    if test $H5_VERS_MINOR -ge 8; then
      h5_use_16=true
    fi
  fi
  if test -z "$new_H5Sselect_hyperslab_ifc"; then
    AC_DEFINE([OLD_H5S_SELECT_HYPERSLAB_IFC], [],
	[whether the old 1.6.3 H5Sselect_hyperslab interface is in use])
  else
    AC_DEFINE([NEW_H5S_SELECT_HYPERSLAB_IFC], [],
	[whether the new H5Sselect_hyperslab interface is in use])
  fi
  if test -n "$h5_use_16"; then
    AC_DEFINE([H5_USE_16_API], [],
	[to force HDF5-1.6 interfaces by default])
  fi
fi
AC_SUBST(H5_VERS_MAJOR)
AC_SUBST(H5_VERS_MINOR)
AC_SUBST(H5_VERS_RELEASE)
HDF5_VERSION=$H5_VERS_MAJOR.$H5_VERS_MINOR.$H5_VERS_RELEASE
AC_SUBST(HDF5_VERSION)

dnl ######################################################################
dnl
dnl  h5diff
dnl
dnl ######################################################################

if test $ac_cv_have_hdf5 = yes; then
  HDF5_BINPATH=`dirname $HDF5_LIBDIR`/bin
  AC_PATH_PROGS(H5DIFF, h5diff, "", $HDF5_BINPATH)
  if test -z "$H5DIFF"; then
    if test -n "$MUST_HAVE_H5DIFF"; then
      AC_MSG_ERROR(h5diff not found in HDF5 tree.  Must reinstall HDF5.)
    fi
  else
    AC_SUBST(H5DIFF)
    if test -z "$new_H5Sselect_hyperslab_ifc"; then
      if test -n "$MUST_HAVE_H5DIFF"; then
        AC_MSG_ERROR(HDF5 must be at version 1.6.4 for regression tests to work.)
      fi
    else
      echo h5diff of version 1.6.4 or greater found.
    fi
  fi
  AC_PATH_PROGS(H5LS, h5ls, "", $HDF5_BINPATH)
  AC_SUBST(H5LS)
  AC_PATH_PROGS(H5DUMP, h5dump, "", $HDF5_BINPATH)
  AC_SUBST(H5DUMP)
else
  if test -n "$MUST_HAVE_H5DIFF"; then
    AC_MSG_ERROR(h5diff not found.  Upgrade hdf5.)
  fi
fi

# use windows path in cygwin so python calls work
case `uname` in
  CYGWIN*)
    H5DIFF=`cygpath -wma "$H5DIFF"`
    H5LS=`cygpath -wma "$H5LS"`
    H5DUMP=`cygpath -wma "$H5DUMP"`
    ;;
esac

dnl ######################################################################
dnl
dnl Determine whether szip needed.
dnl JRC 24jul09: One should not do this here!  Add tx_sz.m4 to your configure.
dnl
dnl ######################################################################

HDF5_NEED_SZIP=`grep H5_HAVE_FILTER_SZIP $HDF5_INCDIR/H5pubconf.h | grep define`

dnl ######################################################################
dnl
dnl    Write to summary file if defined
dnl
dnl ######################################################################

if test -n "$config_summary_file"; then
  if test $ac_cv_have_hdf5 = yes; then
    TX_PRINT_VAR(H5DIFF)
    TX_PRINT_VAR(H5LS)
    TX_PRINT_VAR(H5DUMP)
    if test "$ac_use_fortran_hdf5" = yes; then
      echo USING HDF5 FORTRAN LIBS >> $config_summary_file
    fi
    TX_PRINT_VAR(HDF5_VERSION)
  fi
fi

# For a parallel build, we'd still like vpdatanal built serially. The compiler
# is set properly in vpdatanal/Makefile.am, but even with a serial compiler
# it will try to use parallel HDF5. The vars below are used by
# vpdatanal/Makefile.am to force the use of serial HDF5. If it doesn't exist,
# the parallel version will be used, as before.

HDF5_INCDIR_SERIAL=`echo $HDF5_INCDIR | sed 's/mpi//g'`
HDF5_LIBDIR_SERIAL=`echo $HDF5_LIBDIR | sed 's/mpi//g'`
HDF5_LIBS_SERIAL=`echo $HDF5_LIBS | sed 's/mpi//g'`

if test ! -d "$HDF5_LIBDIR_SERIAL" -o ! -d "$HDF5_INCDIR_SERIAL"; then
  cat <<END
Failed to find serial build of HDF5 in obvious location. Setting _SERIAL vars
to be the same as the parallel vars so as not to break the vpdatanal builds.
END
  HDF5_INCDIR_SERIAL=$HDF5_INCDIR
  HDF5_LIBDIR_SERIAL=$HDF5_LIBDIR
  HDF5_LIBS_SERIAL=$HDF5_LIBS
fi

dnl JRC: no need to do these again
AC_SUBST(HDF5_INCDIR_SERIAL)
AC_SUBST(HDF5_LIBDIR_SERIAL)
AC_SUBST(HDF5_LIBS_SERIAL)

echo "Leaving tx_hdf5.m4."

