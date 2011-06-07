dnl ######################################################################
dnl
dnl File:	cfitsio.m4
dnl
dnl Purpose:	Determine where CCfits is.
dnl
dnl Version:	$Id: cfitsio.m4 3366 2010-01-15 18:43:10Z dws $
dnl
dnl Copyright 2007-2010, Tech-X Corporation. Redistribution allowed provided
dnl this copyright statement remains intact.
dnl
dnl ######################################################################

AC_ARG_WITH(CFITSIO_HOME,
[  --with-CFITSIO_HOME=<cfitsio-home-dir>      to set location of cfitsio],
CFITSIO_HOME="$withval")

dnl
dnl If not known, check in typical directories
dnl
if test -n "$CFITSIO_HOME"; then
  echo "CFITSIO_HOME set by user"
else
  CFITSIO_HOME=/usr/local/cfitsio
fi

CFITSIO_INCDIR=$CFITSIO_HOME/include
CFITSIO_LIBDIR=$CFITSIO_HOME/lib
CFITSIO_LIBS="-L$CFITSIO_LIBDIR -lcfitsio"

if test -z "$CFITSIO_HOME"; then
  AC_MSG_WARN(Unable to find cfitsio.)
else
  AC_DEFINE(HAVE_CFITSIO, , Define if cfitsio is there)
fi

AC_SUBST(CFITSIO_HOME)
AC_SUBST(CFITSIO_INCDIR)
AC_SUBST(CFITSIO_LIBDIR)
AC_SUBST(CFITSIO_LIBS)
