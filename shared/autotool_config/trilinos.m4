dnl ######################################################################
dnl
dnl File:	trilinos.m4
dnl
dnl Purpose:	Determine the locations of trilinos includes and libraries.
dnl
dnl
dnl Version: $Id: trilinos.m4 3517 2010-05-05 18:20:49Z cary $
dnl
dnl Tech-X configure system
dnl
dnl ######################################################################

dnl Include functions from txsearch.m4
builtin(include, config/txsearch.m4)

dnl ######################################################################
dnl
dnl find the location of the trilinos library
dnl
dnl ######################################################################

for i in `echo $SUPRA_SEARCH_PATH | tr ':' ' '`; do
  TRILINOS_SERPATH="$TRILINOS_SERPATH:$i/trilinos"
done
for i in `echo $SUPRA_SEARCH_PATH | tr ':' ' '`; do
  TRILINOS_PARPATH="$TRILINOS_PARPATH:$i/trilinos-par:$i/trilinosmpi"
done

# jrc 2apr10: no longer add default paths here
if test "$parallel" = yes; then
  TRILINOS_PATH=$TRILINOS_PARPATH
#:$abs_top_builddir/../txmodules/trilinosmpi:$HOME/$UNIXFLAVOR/trilinosmpi:/usr/local/trilinosmpi
else
  TRILINOS_PATH=$TRILINOS_SERPATH
#:$abs_top_builddir/../txmodules/trilinos:$HOME/$UNIXFLAVOR/trilinos:/usr/local/trilinos
fi

TX_LOCATE_PKG(
  [trilinos],
  [$TRILINOS_PATH],
  [az_aztec.h],
  [aztecoo, ml, zoltan, amesos, ifpack, epetraext, galeri, triutils, epetra, teuchos],
  [],
  [lib])

if test -n "$TRILINOS_INC_AZ_AZTEC_H"; then
  ac_cv_have_trilinos=yes
else
  ac_cv_have_trilinos=no
fi
if test -n "$TRILINOS_LIB_AMESOS"; then
  AC_DEFINE([HAVE_AMESOS], [], ["Whether Amesos was found"])
fi
if test -n "$TRILINOS_LIB_EPETRAEXT"; then
  AC_DEFINE([HAVE_EPETRAEXT], [], ["Whether EpetraExt was found"])
fi

dnl ######################################################################
dnl
dnl If libraries found, reset variables.  Set version
dnl
dnl ######################################################################

AM_CONDITIONAL(HAVE_TRILINOS, test -n "$TRILINOS_LIBS")
if test -n "$TRILINOS_LIBS"; then
  AC_DEFINE(HAVE_TRILINOS, [], [Defined if trilinos was found.])

  if test -f $TRILINOS_INCDIR/Trilinos_version.h; then
    AC_MSG_CHECKING(trilinos version in $TRILINOS_INCDIR/Trilinos_version.h)
    TRILINOS_MAJOR_VERSION=`grep TRILINOS_MAJOR_VERSION $TRILINOS_INCDIR/Trilinos_version.h | sed -e 's/^.*TRILINOS_MAJOR_VERSION *//' -e 's/ *$//'`
    AC_MSG_RESULT($TRILINOS_MAJOR_VERSION)
    case $TRILINOS_MAJOR_VERSION in
      10.*)
        AC_DEFINE(HAVE_TRILINOS10, [], [Defined if trilinos is version 10.])
        ;;
    esac
  fi

fi
AM_CONDITIONAL(HAVE_TRILINOS10, test "$TRILINOS_MAJOR_VERSION" = 10)

dnl ######################################################################
dnl
dnl Define for whether trilinos was found.
dnl If so, determine its build variables and replace the fortran vars.
dnl
dnl ######################################################################

if test $ac_cv_have_trilinos = yes; then
  if test -f $TRILINOS_INCDIR/Makefile.export.aztecoo.macros; then
    AZTECOO_F77=`grep aztecoo_F77 $TRILINOS_INCDIR/Makefile.export.aztecoo.macros | sed 's/^.* = //'`
    if test -n "$AZTECOO_F77"; then
      AZTECOO_FLIBS=`grep aztecoo_FLIBS $TRILINOS_INCDIR/Makefile.export.aztecoo.macros | sed 's/^.* = //'`
      AZTECOO_LDFLAGS=`grep aztecoo_LDFLAGS $TRILINOS_INCDIR/Makefile.export.aztecoo.macros | sed 's/^.* = //'`
      AZTECOO_LAPACK_LIBS=`grep aztecoo_LAPACK_LIBS $TRILINOS_INCDIR/Makefile.export.aztecoo.macros | sed 's/^.* = //'`
      AZTECOO_BLAS_LIBS=`grep aztecoo_BLAS_LIBS $TRILINOS_INCDIR/Makefile.export.aztecoo.macros | sed 's/^.* = //'`
    else
      AZTECOO_F77=`grep AZTECOO_F77 $TRILINOS_INCDIR/Makefile.export.aztecoo.macros | sed 's/^.* = //'`
      AZTECOO_FLIBS=`grep AZTECOO_FLIBS $TRILINOS_INCDIR/Makefile.export.aztecoo.macros | sed 's/^.* = //'`
      AZTECOO_LDFLAGS=`grep AZTECOO_LDFLAGS $TRILINOS_INCDIR/Makefile.export.aztecoo.macros | sed 's/^.* = //'`
      AZTECOO_LAPACK_LIBS=`grep AZTECOO_LAPACK_LIBS $TRILINOS_INCDIR/Makefile.export.aztecoo.macros | sed 's/^.* = //'`
      AZTECOO_BLAS_LIBS=`grep AZTECOO_BLAS_LIBS $TRILINOS_INCDIR/Makefile.export.aztecoo.macros | sed 's/^.* = //'`
    fi
    AZTECOO_LAPACKBLAS_LIBS="$AZTECOO_LAPACK_LIBS $AZTECOO_BLAS_LIBS"
  elif test -f $TRILINOS_INCDIR/Makefile.export.AztecOO; then
    AZTECOO_LAPACKBLAS_LIBS=`grep AZTECOO_TPL_LIBRARIES $TRILINOS_INCDIR/Makefile.export.AztecOO | sed 's/^.*= //'`
  else
    echo Extra Trilinos libraries not found.
  fi
fi

dnl If we're linking statically, all of the above libs needs to be .a rather
dnl than .so.*. Also removing -lgcc_s if building statically.

if test "$enable_staticlink" = yes; then
    AZTECOO_FLIBS=`echo "$AZTECOO_FLIBS" | \
		sed -e 's/\.so.[[^ ]]*/.a/g' -e 's/-lgcc_s//'`
    AZTECOO_LAPACK_LIBS=`echo "$AZTECOO_LAPACK_LIBS" | \
		sed -e 's/\.so\.[[^ ]]*/.a/g' -e 's/-lgcc_s//'`
    AZTECOO_BLAS_LIBS=`echo "$AZTECOO_BLAS_LIBS" | \
		sed -e 's/\.so\.[[^ ]]*/.a/g' -e 's/-lgcc_s//'`
    AZTECOO_LAPACKBLAS_LIBS=`echo "$AZTECOO_LAPACKBLAS_LIBS" | \
		sed -e 's/\.so\.[[^ ]]*/.a/g' -e 's/-lgcc_s//'`
    AZTECOO_LDFLAGS=`echo "$AZTECOO_LDFLAGS" | sed 's/-lgcc_s//'`
fi

dnl remove extra dirs for aix
case $host in

  *-*-aix*)
    for i in $AZTECOO_FLIBS; do
      isdir=`echo $i | grep -- '^-L'`
      if test -z "$isdir"; then
        newflibs="$newflibs $i"
      fi
    done
    AZTECOO_FLIBS=$newflibs
    ;;

esac
# Whether AZTECOO_FLIBS contains anything
AM_CONDITIONAL(HAVE_AZTECOO_FLIBS, test -n "$AZTECOO_FLIBS")
# havelapackblas=`echo $AZTECOO_LAPACKBLAS_LIBS | egrep -q -- '(^| )-llapack($| )'`
havelapackblas=`echo $AZTECOO_LAPACKBLAS_LIBS | egrep -q -- ' -llapack '`
AM_CONDITIONAL(HAVE_AZTECOO_LAPACKBLAS_LIBS, test -n "$havelapackblas")

dnl ######################################################################
dnl
dnl Now can compose the lapack libraries
dnl
dnl ######################################################################

TRILINOS_LIBS="$TRILINOS_LIBS $AZTECOO_LAPACKBLAS_LIBS $AZTECOO_FLIBS $AZTECOO_LDFLAGS"
# Must clean again
TX_CLEAN_LIBS([TRILINOS_LIBS])

dnl ######################################################################
dnl
dnl Print out modifications
dnl
dnl ######################################################################

if test -n "$config_summary_file"; then
  if test $ac_cv_have_trilinos = yes; then
    echo "  "MODIFICATIONS:  >> "$config_summary_file"
    TX_PRINT_VAR(TRILINOS_LIBS)
    TX_PRINT_VAR(TRILINOS_RPLIBS)
    TX_PRINT_VAR(TRILINOS_LTLIBS)
    TX_PRINT_VAR(TRILINOS_ALIBS)
    TX_PRINT_VAR(TRILINOS_MAJOR_VERSION)
  else
    echo "Trilinos not found." >> $config_summary_file
  fi
fi

AC_SUBST(TRILINOS_LIBS)

