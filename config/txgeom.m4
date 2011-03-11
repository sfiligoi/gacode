dnl ######################################################################
dnl
dnl File:      	txgeom.m4
dnl
dnl Purpose:	Determine where the txgeom files are.
dnl
dnl Version:	$Id: txgeom.m4 3725 2010-10-27 17:16:02Z mmiah $
dnl
dnl Copyright 2001-2010, Tech-X Corporation.  Redistribution allowed provided
dnl this copyright statement remains intact.
dnl
dnl ######################################################################

dnl ######################################################################
dnl
dnl Determine the location of txgeom
dnl
dnl ######################################################################

builtin(include, config/txsearch.m4)

dnl Locate all the txgeom related includes and libraries
AC_ARG_WITH(txgeom-subdir,
  AC_HELP_STRING([--with-txgeom-subdir=<txgeom-subdir>],
    [to look for installations under <txgeom-subdir>]),
    TXGEOM_SUBDIR="$withval")
if test -n "$TXGEOM_SUBDIR"; then
  TXGEOM_SUBDIRS=$TXGEOM_SUBDIR
else
  TXGEOM_SUBDIRS=txgeom
fi

unset TXGEOM_PATH

for i in `echo $SUPRA_SEARCH_PATH | tr ':' ' '`; do
  for j in $TXGEOM_SUBDIRS; do
    TXGEOM_PATH="$TXGEOM_PATH:$i/$j"
  done
done

TX_LOCATE_PKG(
  [txgeom],
  [$TXGEOM_PATH],
  [TxGeomLibrary.h],
  [txgeomlibrary, txgeommesh, txgeomgeometry, txgeomtopology, txgeombase],
  [],
  [])
