dnl ######################################################################
dnl
dnl File:       defaultcomps.m4
dnl
dnl Purpose:    Select the default C++ compiler for a given OS, then the
dnl 		default C and fortran compilers from that.
dnl		E.g. xlc_r for xlC_r.  If parallel is not set, then it
dnl		is set to no.
dnl
dnl Version:    $Id: defaultcomps.m4 3647 2010-09-03 22:20:42Z veitzer $
dnl
dnl Copyright 2001-2010, Tech-X Corporation.
dnl Redistribution allowed provided this copyright statement remains intact.
dnl
dnl ######################################################################

AC_DEFUN([TX_DEFAULT_CXX], [
  if test -z "$CXX"; then
    if test -z "$parallel"; then
      parallel=no
    fi
    case $host in
      *-*-aix*)
        case $parallel in
          no)  CXX=xlC_r;;
          yes) CXX=mpCC_r;;
        esac
        ;;
      *)
        case $parallel in
          no)  CXX=g++;;
          yes) CXX=mpicxx;;
        esac
        ;;
    esac
  fi
])


AC_DEFUN([TX_DEFAULT_CC], [
  if test -z "$CC"; then
    if test -n "$CXX"; then
      CXXBASE=`basename $CXX`
      case "$CXXBASE" in
        g++ | c++) CC=gcc ;;
        mpiCC | mpicxx) CC=mpicc ;;
        mpCC_r) CC=mpcc_r ;;
        aCC) CC=cc ;;
        CC) CC=cc ;;
        xlC) CC=xlc ;;
        xlC_r) CC=xlc_r ;;
        *) CC=gcc ;;
      esac
      ABSCXX=`which $CXX`
      CXXDIR=`dirname $ABSCXX`
      if test -x $CXXDIR/$CC; then
        CC=$CXXDIR/$CC
      fi
    else
      CC=gcc
    fi
  fi
])

AC_DEFUN([TX_DEFAULT_FC], [
  if test -z "$FC"; then
    CCBASE=`basename $CC`
    case "$CCBASE" in
      gcc) FC=gfortran ;;
      mpicc) FC=mpif90 ;;
      mpcc) FC=mpxlf ;;
      mpcc_r) FC=mpxlf_r ;;
      acc) FC=f90 ;;
      cc) FC=ftn ;;
      xlc) FC=xlf ;;
      xlc_r) FC=xlf_r ;;
      *) FC=f90 ;;
    esac
    ABSCC=`which $CC`
    CCDIR=`dirname $ABSCC`
# Make name absolute and consistent with CC if possible
    if test -x $CCDIR/$FC; then
      FC=$CCDIR/$FC
    fi
    if ! which $FC 1>/dev/null 2>&1; then
      unset FC
    fi
  fi
])


# Default fixed form Fortran
# Usually we are just setting F77 for fixed format code
# (automake rules), so just use the FC compiler
# Anything else becomes very dangerous in practice.
AC_DEFUN([TX_DEFAULT_F77], [
  if test -z "$F77"; then
   F77=$FC
   AC_SUBST(F77)
  fi
])


