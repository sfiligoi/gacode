dnl ######################################################################
dnl
dnl File:       txsuperlu.m4
dnl
dnl Purpose:    Determine if SuperLU is installed.
dnl
dnl Version:    $Id: tx_superlu.m4 3710 2010-10-13 19:52:43Z cary $
dnl
dnl Copyright 2001-2010, Tech-X Corporation.  Redistribution allowed
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
  SUPERLU_SP=$HOME:$HOME/$UNIXFLAVOR:$HOME/software:$SUPRA_SEARCH_PATH
  unset SUPERLU_SEQ_PATH
  for i in `echo $SUPERLU_SP | tr ':' ' '`; do
    SUPERLU_SEQ_PATH="$SUPERLU_SEQ_PATH:$i/superlu_seq:$i/petsc"
  done
  SUPERLU_SEQ_PATH="$SUPERLU_SEQ_PATH:/usr/common/acts/SuperLU/SuperLU_3.0"
  echo SUPERLU_SEQ_PATH = $SUPERLU_SEQ_PATH
  if test -n "$PETSC_INC_PETSC_H"; then
    PETSC_INCDIR1=`dirname $PETSC_INC_PETSC_H`
    PETSC_DIR=`dirname $PETSC_INCDIR1`
    SUPERLU_SEQ_PATH=$PETSC_DIR:$SUPERLU_SEQ_PATH
  fi
  TX_LOCATE_PKG(
	[SUPERLU_SEQ],
	[$SUPERLU_SEQ_PATH],
	[colamd.h,old_colamd.h,slu_Cnames.h,slu_ddefs.h,slu_util.h],
	[superlu_4.0],
	[include:SRC],
	[lib:.],)

  if test -n "$SUPERLU_SEQ_PATH_LIBDIR"; then
    symbol="superlu_free"
    AC_CHECK_LIB(
	[superlu_4.0],
	[$symbol],,
	AC_MSG_ERROR($symbol not found in libsuperlu_4.0.),
    )
  fi
fi

AM_CONDITIONAL([HAVE_SUPERLU_SEQ], [test "$FOUND_SUPERLU_SEQ" != "no"])

dnl ##################################################################
dnl
dnl Locate the distributed package, verify it, set conditional
dnl
dnl ##################################################################

if test "$WANT_SUPERLU_DIST" = yes; then

  SUPERLU_DIST_PATH=$abs_top_builddir/../txmodules/superlumpi:$HOME/superlu_dist:$HOME/$UNIXFLAVOR/superlumpi:$HOME/software/superlu-mpi:/usr/local/superlu_dist:/usr/local/superlumpi
  TX_LOCATE_PKG(
	[SUPERLU_DIST],
	[$SUPERLU_DIST_PATH],
	[colamd.h,old_colamd.h,slu_Cnames.h,slu_ddefs.h,slu_util.h,supermatrix.h],
	[superlu_4.0])

  if test -n "$SUPERLU_DIST_PATH_LIBDIR"; then
    symbol="superlu_free"
    AC_CHECK_LIB(
	[superlu_4.0],
	[$symbol],,
	AC_MSG_ERROR($symbol not found in libsuperlu_4.0.),
    )
  fi
fi

AM_CONDITIONAL([HAVE_SUPERLU_DIST], [test "$FOUND_SUPERLU_DIST" != "no"])


