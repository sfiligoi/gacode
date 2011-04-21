dnl ######################################################################
dnl
dnl File:	tx_sz.m4
dnl
dnl Purpose:	Determine the locations of sz includes and libraries
dnl
dnl Version: $Id: tx_sz.m4 3450 2010-03-27 21:11:39Z cary $
dnl
dnl Tech-X configure system
dnl
dnl This is being split off from hdf5.m4, so it can be called
dnl after RPATH_FLAGS is defined.
dnl
dnl    SZIP_INCDIR
dnl    SZIP_INC
dnl    SZIP_LIBDIR
dnl    SZIP_LIBS
dnl
dnl ######################################################################

dnl Includes functions from txsearch.m4
builtin(include, shared/autotool_config/txsearch.m4)

dnl ######################################################################
dnl
dnl Construct search path for szip
dnl
dnl ######################################################################

# Look first for szip consistent with hdf5.
if test -n "$HDF5_LIBDIR"; then
  SZ_PATH="$SZ_PATH:$HDF5_LIBDIR"
fi
# Add directories from supra search path
if test -n "$SUPRA_SEARCH_PATH"; then
  for i in `echo $SUPRA_SEARCH_PATH | tr ':' ' '`; do
    SZ_PATH="$SZ_PATH:$i/szip:$i/hdf5"
  done
fi
# Add directory from module package.  Possible since serial.
if test -n "$SZIP_DIR"; then
  SZ_PATH="$SZ_PATH:$SZIP_DIR"
fi
# Add other default directories
SZ_PATH=$SZ_PATH:$HOME/$UNIXFLAVOR/hdf5:/contrib/hdf5:/usr/local/hdf5

dnl ######################################################################
dnl
dnl Find szip libraries
dnl
dnl ######################################################################

TX_LOCATE_PKG(
    [SZ],
    [$SZ_PATH],
    [szlib.h],
    [sz])
TX_CLEAN_LIBS([SZ_LIBS])
dnl AC_SUBST(SZ_RPLIBS)
dnl AC_SUBST(SZ_LTLIBS)
dnl AC_SUBST(SZ_ALIBS)

# Needed because set blank if not found and still okay
# echo SZ_LIBS = $SZ_LIBS
if test -n "$SZ_LIBDIR"; then
  SZ_LIB=-lsz
  SZ_LDFLAG=-L$SZ_LIBDIR
  SZ_RLDFLAG=${RPATH_FLAG}${SZ_LIBDIR}
fi
AC_SUBST(SZ_LIB)
AC_SUBST(SZ_LDFLAG)
AC_SUBST(SZ_RLDFLAG)

dnl ######################################################################
dnl
dnl Fix the RPATH_FLAG is possible
dnl
dnl ######################################################################

dnl mmiah: why is this even necessary?  TX_LOCATE_PKG should be adding
dnl        the rpath flags, and if not, it needs to be fixed.  Leaving
dnl        this here for now, as I do not have an SZ testbed to confirm
dnl JRC: legacy.  It used to be that szlib was found at the very
dnl      beginning to determine parallelism, when RPATH_FLAG was set to
dnl 	 $(RPATH_FLAG) to be determined later.  However, in most cases
dnl      now, tx_sz.m4 is called later, after RPATH_FLAG is set, and so
dnl      the value can be inserted, making a smaller burden on developers.
dnl      Of course, this logic could be added to TX_LOCATE_PKG.
if test -z "$RPATH_FLAG"; then
  AC_MSG_WARN(RPATH_FLAG not defined.  Move tx_sz.m4 after libs.m4.)
else
  SZ_LIBS=`echo $SZ_LIBS | sed -e 's/\$(//'  -e 's/)//' -e "s/RPATH_FLAG/$RPATH_FLAG/"`
  SZ_RPLIBS=`echo $SZ_RPLIBS | sed -e 's/\$(//'  -e 's/)//' -e "s/RPATH_FLAG/$RPATH_FLAG/"`
fi
if test -n "$config_summary_file" -a x"$SZ_DOSEARCH" != xno; then
  echo "  SZ mods, additional:"            >> $config_summary_file
  # echo "    SZ_LIB:     $SZ_LIB"       >> $config_summary_file
  # echo "    SZ_LIBS:    $SZ_LIBS"      >> $config_summary_file
  # echo "    SZ_LDFLAG:  $SZ_LDFLAG"    >> $config_summary_file
  # echo "    SZ_RLDFLAG: $SZ_RLDFLAG"   >> $config_summary_file
  TX_PRINT_VAR(SZ_LIB)
  TX_PRINT_VAR(SZ_LIBS)
  TX_PRINT_VAR(SZ_RPLIBS)
  dnl TX_PRINT_VAR(SZ_LTLIBS)
  dnl TX_PRINT_VAR(SZ_ALIBS)
  TX_PRINT_VAR(SZ_LDFLAG)
  TX_PRINT_VAR(SZ_RDFLAG)
fi

