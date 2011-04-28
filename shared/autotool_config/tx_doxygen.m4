dnl ######################################################################
dnl
dnl File:	tx_doxygen.m4
dnl
dnl Purpose:	Determine the location of doxygen and set the directory
dnl		for its output
dnl
dnl Version: $Id: tx_doxygen.m4 3715 2010-10-15 20:23:39Z kruger $
dnl
dnl ######################################################################

dnl ######################################################################
dnl
dnl Look for the executable.
dnl
dnl ######################################################################

AC_PATH_PROGS(DOXYGEN, doxygen, "")
if test -z "$DOXYGEN"; then
  AC_MSG_WARN(The doxygen utility was not found in your path.  Combined C++ documentation will not be made.)
fi
AC_SUBST(DOXYGEN)

dnl ######################################################################
dnl
dnl Allow setting of doxygen dir.  Otherwise defaults to docs/cxxapi
dnl
dnl ######################################################################

dnl  We are now using doxygen for f90 codes as well so we need the
dnl  configure.ac to specify a different default at times
if test -z "$ac_default_doxygen_dir"; then
   ac_default_doxygen_dir=docs/cxxapi
fi

AC_ARG_WITH(doxygen-docsdir,
	[  --with-doxygen-docsdir=<where to put doxygen output>],
	DOXYGEN_DOCSDIR="$withval", DOXYGEN_DOCSDIR="")
if test -z "$DOXYGEN_DOCSDIR"; then
  DOXYGEN_DOCSDIR=$ac_default_doxygen_dir
fi
if test ! -d $DOXYGEN_DOCSDIR; then
  mkdir -p $DOXYGEN_DOCSDIR
fi
DOXYGEN_DOCSDIR=`(cd $DOXYGEN_DOCSDIR; pwd)`
AC_SUBST(DOXYGEN_DOCSDIR)
DOXYGEN_SUBDIR=`basename $DOXYGEN_DOCSDIR`
AC_SUBST(DOXYGEN_SUBDIR)
if test -n "$config_summary_file"; then
  echo                                            >> $config_summary_file
  echo "Documentation built with" >> $config_summary_file
  # echo "  DOXYGEN:  $DOXYGEN"  >> $config_summary_file
  TX_PRINT_VAR(DOXYGEN)
  TX_PRINT_VAR(DOXYGEN_DOCSDIR)
fi

