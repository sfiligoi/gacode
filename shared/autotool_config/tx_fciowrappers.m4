dnl ######################################################################
dnl
dnl File:	tx_fciowrappers.m4
dnl
dnl Purpose:	Find the fciowrappers package and grab the location
dnl of the HDF5 and NETCDF directories
dnl
dnl Version: $Id: tx_fciowrappers.m4 3806 2011-03-07 21:31:56Z kruger $
dnl
dnl Tech-X configure system
dnl
dnl ######################################################################

dnl ######################################################################
dnl
dnl Check --with setting
dnl
dnl ######################################################################

AC_ARG_WITH(fciowrappers-dir, AC_HELP_STRING([--with-fciowrappers-dir],
        [fciowrappers installation directory]),
        FCIOWRAPPERS_DIR="$withval")

dnl ######################################################################
dnl
dnl Set the search path.
dnl
dnl ######################################################################

FCIO_SP=$SUPRA_SEARCH_PATH
FCIO_SERPATH=$FCIOWRAPPERS_DIR
for i in `echo $FCIO_SP | tr ':' ' '`; do
  FCIO_SERPATH="$FCIO_SERPATH:$i/fciowrappers"
done
FCIO_PARPATH=$FCIOWRAPPERS_DIR
for i in `echo $FCIO_SP | tr ':' ' '`; do
  FCIO_PARPATH="$FCIO_PARPATH:$i/fciowrappers-par:$i/fciowrappersmpi"
done

if test -z "$PAR_FCIO_IS_DEFAULT"; then
  if test -n "$back_end_node"; then
    PAR_FCIO_IS_DEFAULT=yes
  elif test -n "$parallel"; then
    PAR_FCIO_IS_DEFAULT=$parallel
  elif test -n "$paralleldefault"; then
    PAR_FCIO_IS_DEFAULT=$paralleldefault
  else
    PAR_FCIO_IS_DEFAULT=yes
  fi
fi
if test "$PAR_FCIO_IS_DEFAULT" = yes; then
  FCIO_PATH=$FCIO_PATH:$FCIO_PARPATH:$FCIO_SERPATH
else
  FCIO_PATH=$FCIO_PATH:$FCIO_SERPATH:$FCIO_PARPATH
fi

dnl ######################################################################
dnl
dnl Search.
dnl
dnl ######################################################################

if test "$USE_FCIO_ONLY" = netcdf; then
TX_LOCATE_PKG(
  [FCIOWRAPPERS],
  [$FCIO_PATH],
  [vshdf5_dummy.h],
  [ezcdf],
  [include],
  [lib/$FC_LIBSUBDIR],
)
elif test "$USE_FCIO_ONLY" = hdf5; then
TX_LOCATE_PKG(
  [FCIOWRAPPERS],
  [$FCIO_PATH],
  [vshdf5_dummy.h],
  [vshdf5],
  [include],
  [lib/$FC_LIBSUBDIR],
)
else
TX_LOCATE_PKG(
  [FCIOWRAPPERS],
  [$FCIO_PATH],
  [vshdf5_dummy.h],
  [ezcdf, vshdf5],
  [include],
  [lib/$FC_LIBSUBDIR],
)
fi

dnl ######################################################################
dnl
dnl  TX_LOCATE_PKG is an all or nothing about whether it is found
dnl
dnl ######################################################################

if test -n "$FCIOWRAPPERS_LIBS"; then
  FOUND_FCIOWRAPPERS=yes
fi

dnl ######################################################################
dnl
dnl Grab the HDF5 and netCDF information if requested and found.
dnl If requested, not found, unset AM_CONDITIONALs, HAVE_HDF5 and HAVE_NETCDF.
dnl If not requested, assume previously found and check consistency.
dnl
dnl ######################################################################
dnl     FOUND_pkgname:      set to "yes" if package was found, or "no"

IIICS=$FCIOWRAPPERS_INCDIR/../share/config.summary
if test "$USE_IO_INFO" = true; then
  if test "$FOUND_FCIOWRAPPERS" == "yes"; then
    FCIOWRAPPERS_VARS="HDF5_INC HDF5_INCDIR HDF5_LIBDIR HDF5_LIBS HDF5_RPLIBS HDF5_LTLIBS HDF5_ALIBS HDF5_MOD HDF5_MODDIR HDF5_FC_INC HDF5_FC H5DIFF H5DUMP H5LS HDF5_VERSION GPFS_LIB GPFS_LIBDIR GPFS_LIBS NETCDF_INC NETCDF_INCDIR NETCDF_LIBDIR NETCDF_LIBS NETCDF_ALIBS NETCDF_RPLIBS SZ_INC SZ_INCDIR SZ_INC_SZLIB_H SZ_LDFLAG SZ_LIB_SZ SZ_LIBDIR SZ_LIBS SZ_RPLIBS SZ_LTLIBS SZ_ALIBS"
    for i in $FCIOWRAPPERS_VARS; do
      # line=`grep ' '$i' ' $IIICS`
      value=`grep ' '$i' ' $IIICS | sed -e 's/^.*: *//'`
      # echo Value for $i is $value.
# Use single quotes for $(RPATH_FLAG)
      eval $i=\'$value\'
    done
    AM_CONDITIONAL(HAVE_HDF5, test "$HDF5_LIBS" != "" )
    AM_CONDITIONAL(HAVE_NETCDF, test "$NETCDF_LIBS" != "" )
    if test -n "$NETCDF_LIBS"; then
      AC_DEFINE(HAVE_NETCDF, [], [Defined if netcdf found by fciowrappers])
    fi
    if test -n "$HDF5_LIBS"; then
      AC_DEFINE(HAVE_HDF5, [], [Defined if hdf5 found by fciowrappers])
    fi
    if test "$HDF5_LIBS" != ""; then
       case "$SERIALFC" in
	 *xlf*) FPPFLAGS="-WF,-D__HAVE_VSHDF5 $FPPFLAGS";;
         *)    FPPFLAGS="-D__HAVE_VSHDF5 $FPPFLAGS";;
       esac
       echo "   Using vshdf5. "                                      >> $config_summary_file
       echo "   Type: 'grep HDF5 Makefile' to see HDF5 location"     >> $config_summary_file
    fi
    if test "$NETCDF_LIBS" != ""; then
       case "$SERIALFC" in
	 *xlf*) FPPFLAGS="-WF,-D__HAVE_EZCDF $FPPFLAGS";;
         *)    FPPFLAGS="-D__HAVE_EZCDF $FPPFLAGS";;
       esac
       echo "   Using ezcdf. "                                       >> $config_summary_file
       echo "   Type: 'grep NETCDF Makefile' to see NETCDF location" >> $config_summary_file
       echo                                                          >> $config_summary_file
    fi
  else
    AM_CONDITIONAL(HAVE_HDF5, test 1 == 0 )
    AM_CONDITIONAL(HAVE_NETCDF, test 1 == 0 )
    dnl SEK: I'm not sure whether I really need this.  I am not an
    dnl AC_DEFINE expert, but I don't think this hurts.
    AC_DEFINE(HAVE_NETCDF, 0, [Defined if netcdf found by fciowrappers])
    AC_DEFINE(HAVE_HDF5, 0, [Defined if hdf5 found by fciowrappers])
  fi
  for i in $FCIOWRAPPERS_VARS; do
    TX_PRINT_VAR($i)
  done
  AC_SUBST(HDF5_INC)
  AC_SUBST(HDF5_INCDIR)
  AC_SUBST(HDF5_LIBDIR)
  AC_SUBST(HDF5_LIBS)
  AC_SUBST(HDF5_RPLIBS)
  AC_SUBST(HDF5_LTLIBS)
  AC_SUBST(HDF5_ALIBS)
  AC_SUBST(HDF5_MOD)
  AC_SUBST(HDF5_MODDIR)
  AC_SUBST(HDF5_FC)
  AC_SUBST(HDF5_FC_INC)
  AC_SUBST(H5DIFF)
  AC_SUBST(H5DUMP)
  AC_SUBST(H5LS)
  AC_SUBST(HDF5_VERSION)
  AC_SUBST(GPFS_LIB)
  AC_SUBST(GPFS_LIBDIR)
  AC_SUBST(GPFS_LIBS)
  AC_SUBST(NETCDF_INC)
  AC_SUBST(NETCDF_INCDIR)
  AC_SUBST(NETCDF_LIBDIR)
  AC_SUBST(NETCDF_LIBS)
  AC_SUBST(NETCDF_ALIBS)
  AC_SUBST(NETCDF_RPLIBS)
  AC_SUBST(SZ_INC)
  AC_SUBST(SZ_INCDIR)
  AC_SUBST(SZ_INC_SZLIB_H)
  AC_SUBST(SZ_LDFLAG)
  AC_SUBST(SZ_LIB_SZ)
  AC_SUBST(SZ_LIBDIR)
  AC_SUBST(SZ_LIBS)
  AC_SUBST(SZ_RPLIBS)
  AC_SUBST(SZ_LTLIBS)
  AC_SUBST(SZ_ALIBS)
else
  if test -n "$HDF5_LIBDIR"; then
    FCIO_HDF5_LIBDIR=`grep ' HDF5_LIBDIR ' $IIICS | sed -e 's/^.*: *//'`
    FCIO_HDF5_LIBDIR=`(cd $FCIO_HDF5_LIBDIR; pwd -P)`
    if test "$FCIO_HDF5_LIBDIR" = "$HDF5_LIBDIR"; then
      echo fciowrappers uses consistent hdf5.
    else
      AC_MSG_ERROR(fciowrappers inconsistent with hdf5.  HDF5_LIBDIR = $HDF5_LIBDIR.  FCIO_HDF5_LIBDIR = $FCIO_HDF5_LIBDIR)
    fi
  else
    AC_MSG_ERROR(fciowrappers found but not hdf5.  Place latter first.)
  fi
fi

