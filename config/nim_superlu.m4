dnl ######################################################################
dnl
dnl File:       txsuperlu.m4
dnl
dnl Purpose:    Determine if SuperLU is installed.
dnl
dnl Version:    $Id: nim_superlu.m4 3464 2010-04-05 12:42:26Z cary $
dnl
dnl Copyright Tech-X Corporation, 2001-2005.  Redistribution allowed
dnl provided this copyright statement remains intact.
dnl
dnl ######################################################################

dnl ##################################################################
dnl
dnl Include functions from txsearch.m4
dnl
dnl ##################################################################
builtin(include, config/txsearch.m4)

dnl ##################################################################
dnl
dnl Locate the sequential package, verify it, set conditional
dnl
dnl ##################################################################

if test "$WANT_SUPERLU_SEQ" = yes; then

dnl mmiah: the slu_supermatrix.h change is not backwards compatible.
dnl        I am temporarily removing the file from the search list
dnl        until I implement a better solution.

dnl PETSc, then near, then home, then system.
  SUPERLU_SP=$SUPRA_SEARCH_PATH:$HOME/lib
  unset SUPERLU_SEQ_PATH
  for i in `echo $SUPERLU_SP | tr ':' ' '`; do
    SUPERLU_SEQ_PATH="$SUPERLU_SEQ_PATH:$i/superlu"
  done
  SUPERLU_SEQ_PATH="$SUPERLU_SEQ_PATH:/usr/common/acts/SuperLU/SuperLU_2.3"
  SUPERLU_SEQ_PATH="$SUPERLU_SEQ_PATH:/opt:/usr/local"
  SUPERLU_SEQ_PATH="$SUPERLU_SEQ_PATH:$HOME/lib/SuperLU"
  echo SUPERLU_SEQ_PATH = $SUPERLU_SEQ_PATH
  TX_LOCATE_PKG(
	[SUPERLU_SEQ],
	[$SUPERLU_SEQ_PATH],
	[colamd.h,old_colamd.h,slu_Cnames.h,slu_ddefs.h,slu_util.h],
	[superlu],
	[include:SRC],
	[lib:.],)

  if test -n "$SUPERLU_SEQ_PATH_LIBDIR"; then
    symbol="superlu_free"
    AC_CHECK_LIB(
	[superlu],
	[$symbol],,
	AC_MSG_ERROR($symbol not found in libsuperlu.),
    )
  fi
else
AM_CONDITIONAL([HAVE_SUPERLU_SEQ], [test 1 == 0])
fi


dnl ##################################################################
dnl
dnl Locate the distributed package, verify it, set conditional
dnl
dnl ##################################################################

if test "$WANT_SUPERLU_DIST" = yes; then

  SUPERLU_SP=$SUPRA_SEARCH_PATH:$HOME/lib
  unset SUPERLU_DIST_PATH
  for i in `echo $SUPERLU_SP | tr ':' ' '`; do
    SUPERLU_DIST_PATH="$SUPERLU_DIST_PATH:$i/superlu_dist"
  done

  TX_LOCATE_PKG(
	[SUPERLU_DIST],
	[$SUPERLU_DIST_PATH],
	[old_colamd.h,slu_Cnames.h,slu_ddefs.h,slu_util.h,supermatrix.h],
	[superlu,superlu_dist],
	[include:SRC],
	[lib:.],)

  if test -n "$SUPERLU_DIST_PATH_LIBDIR"; then
    symbol="superlu_free"
    AC_CHECK_LIB(
	[superlu],
	[$symbol],,
	AC_MSG_ERROR($symbol not found in libsuperlu.),
    )
  fi
else
AM_CONDITIONAL([HAVE_SUPERLU_DIST], [test 1 == 0])
fi


