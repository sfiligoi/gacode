dnl ######################################################################
dnl
dnl File:	txbase.m4
dnl
dnl Purpose:	Determine the locations of txbase includes and libraries
dnl
dnl Version: 	$Id: txbase.m4 3737 2010-11-03 18:57:29Z cary $
dnl
dnl Tech-X configure system
dnl
dnl ######################################################################

dnl ######################################################################
dnl
dnl Construct path from anything in environment or on command line
dnl followed by supra-search-path directories.  Then locate package.
dnl
dnl ######################################################################

AC_ARG_WITH(txbase-dir, AC_HELP_STRING([--with-txbase-dir],
        [txbase installation directory]),
        TXBASE_DIR="$withval")
TXBASE_PATH=$TXBASE_DIR
TXBASE_SP=$SUPRA_SEARCH_PATH
if test "$parallel" = no; then
  txbase_dirs="txbase"
else
  txbase_dirs="txbase-par txbase-ben txbase"
fi
for j in $txbase_dirs; do
  for i in `echo $TXBASE_SP | tr ':' ' '`; do
    TXBASE_PATH="$TXBASE_PATH:$i/$j"
  done
done
TX_LOCATE_PKG(
  [TXBASE],
  [$TXBASE_PATH],
  [txbase_version.h],
  [txbase, txstd],
  [include],
  [lib/$COMPDIR])
AM_CONDITIONAL(HAVE_TXBASE, test -n "$TXBASE_LIBS")
if test -n "$TXBASE_LIBS"; then
  AC_DEFINE(HAVE_TXBASE, [], [Defined if TXBASE found])
fi

dnl ######################################################################
dnl
dnl Find txbase version file
dnl
dnl ######################################################################

dnl Check TxBase version number
if test -n "$TXBASE_INCDIR"; then
  TX_PATH_FILE(TXBASE_VERSION_H, txbase_version.h, "", $TXBASE_INCDIR)
  if test -z "$TXBASE_VERSION_H"; then
dnl JRC: changed to warning to that configure.ac does the kill.
    AC_MSG_WARN(Cannot find txbase_version.h in $TXBASE_INCDIR!
    Please install TxBase or use --with-txbase-dir=)
  fi
fi

dnl ######################################################################
dnl
dnl Verify that version is sufficiently new
dnl
dnl ######################################################################

if test -n "$TXBASE_VERSION_H"; then
  AC_MSG_CHECKING(TxBase version)
  TXBASE_VERSION=`grep TXBASE_VERSION $TXBASE_VERSION_H | grep -v undef | sed 's/^.* \"//' | sed 's/\"//g'`
  if test -z "$TXBASE_VERSION"; then
    AC_MSG_RESULT(not found.)
    AC_MSG_WARN(TXBASE_VERSION not present in $TXBASE_VERSION_H.  Please reinstall.)
  else
    AC_MSG_RESULT($TXBASE_VERSION)
    TXBASE_MAJOR_VERSION=`echo $TXBASE_VERSION | sed 's/\..*$//'`

    if test -z "$TXBASE_OK_MAJOR_VERSION"; then
      TXBASE_OK_MAJOR_VERSION=2
    fi
    if test -z "$TXBASE_OK_MINOR_VERSION"; then
      TXBASE_OK_MINOR_VERSION=4
    fi
    if test -z "$TXBASE_OK_RELEASE"; then
      TXBASE_OK_RELEASE=0
    fi
    if test -z "$TXBASE_OK_VERSION"; then
      TXBASE_OK_VERSION=$TXBASE_OK_MAJOR_VERSION.$TXBASE_OK_MINOR_VERSION.$TXBASE_OK_RELEASE
    fi

    if test $TXBASE_MAJOR_VERSION -lt $TXBASE_OK_MAJOR_VERSION; then
      AC_MSG_ERROR(Major version must be at least $TXBASE_OK_MAJOR_VERSION)
    fi
    if test $TXBASE_MAJOR_VERSION -eq $TXBASE_OK_MAJOR_VERSION; then
      TXBASE_MINOR_VERSION=`echo $TXBASE_VERSION | sed "s/^$TXBASE_MAJOR_VERSION\.//" | sed 's/\..*$//'`
      # echo TXBASE_MINOR_VERSION = $TXBASE_MINOR_VERSION
      if test $TXBASE_MINOR_VERSION -lt $TXBASE_OK_MINOR_VERSION; then
        AC_MSG_ERROR(Minor version must be at least $TXBASE_OK_MINOR_VERSION)
      fi
      if test $TXBASE_MINOR_VERSION -eq $TXBASE_OK_MINOR_VERSION; then
        TXBASE_RELEASE=`echo $TXBASE_VERSION | sed 's/^.*\.//'`
        # echo TXBASE_RELEASE = $TXBASE_RELEASE
        if test $TXBASE_RELEASE -lt $TXBASE_OK_RELEASE; then
          AC_MSG_ERROR(Minor version must be at least $TXBASE_OK_RELEASE)
        fi
      fi
    fi
    echo TxBase of sufficiently high version
  fi

dnl JRC, 22May06: Why define this?  Have substituted it.
  AC_DEFINE_UNQUOTED([TXBASE_INCDIR], "$TXBASE_INCDIR",
	"txbase include directory")

fi

