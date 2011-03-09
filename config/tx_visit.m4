dnl ######################################################################
dnl
dnl $Id: tx_visit.m4 3573 2010-06-23 04:44:06Z mdurant $
dnl
dnl Find VisIt
dnl
dnl ######################################################################

builtin(include, config/txsearch.m4)

dnl ######################################################################
dnl
dnl Obtain visit hostnaming.  Not complete.  See svn_bin/visit-bin-dist
dnl
dnl ######################################################################

os=`uname -s | tr '[A-Z]' '[a-z]' | tr -d '[0-9]'`
version=`uname -r`
mach=`uname -m`
# Do visit transformation of machine name
case $mach in
  i?86)
    case $os in
      darwin) mach=i386;;
      linux) mach=intel;;
    esac
    ;;
  powerpc) mach=ppc;;
esac
# Default for VISIT_BINDIR followed by special cases
VISIT_OS_MACH=$os-$mach
case $os in
  freebsd) VISIT_OS_MACH=FreeBSD-$version-$mach;;
esac
AC_SUBST(VISIT_OS_MACH)

dnl ######################################################################
dnl
dnl Find directory containing include, lib, etc.
dnl
dnl ######################################################################

AC_ARG_WITH(visit-dir, AC_HELP_STRING([--with-visit-dir=<visit install dir>],
        [set the VisIt installation directory]),
        VISIT_DIR="$withval",)
AC_ARG_WITH(visit-version,
        AC_HELP_STRING([--with-visit-version=<visit version>],
        [set the VisIt version]),
        VISIT_VERSION="$withval", VISIT_VERSION=current)
# Set path
unset VISIT_PATH
if test -n "$VISIT_DIR"; then
  VISIT_PATH=$VISIT_DIR/$VISIT_VERSION/$VISIT_OS_MACH
fi
for i in `echo $SUPRA_SEARCH_PATH | tr ':' ' '`; do
  VISIT_PATH="$VISIT_PATH:$i/visit/$VISIT_VERSION/$VISIT_OS_MACH"
done
echo VISIT_PATH = $VISIT_PATH
TX_LOCATE_PKG(
    [visit],
    [$VISIT_PATH],
    [make-variables, avtSTMDFileFormat.h],
    [visitpy],
    [include:include/visit],
)
if test -z "$VISIT_INCDIR"; then
  AC_MSG_ERROR(VisIt includes not found.  Use --with-visit-dir= to set their location.)
fi

dnl ######################################################################
dnl
dnl Find system plugin directory
dnl
dnl ######################################################################

VISIT_SYSDIR=`dirname $VISIT_LIBDIR`
VISIT_PLUGINDIR=$VISIT_SYSDIR/plugins/databases
AC_SUBST(VISIT_PLUGINDIR)

dnl ######################################################################
dnl
dnl Find executables
dnl
dnl ######################################################################

AC_MSG_CHECKING(for VisIt)
VISIT_TOPINCDIR=`dirname $VISIT_INCDIR`
VISIT_ALLINCDIR=`dirname $VISIT_TOPINCDIR`
vdir=`dirname $VISIT_ALLINCDIR`
VISIT_DIR2=`dirname $vdir`
VISIT_PATH=$VISIT_ALLINCDIR/bin:$VISIT_DIR2/bin
AC_PATH_PROGS(VISIT, visit, "", $VISIT_PATH)
if test -z "$VISIT"; then
  AC_MSG_ERROR(VisIt not found.  Please report this error so it can be fixed.)
fi
VISIT_BINDIR=`dirname $VISIT`
XML2MAKEFILE=$VISIT_BINDIR/xml2makefile
dnl XML2MAKEFILE=$abs_top_builddir/txutils/xml2makefile.sh
XML2MAKEFILE_FLAGS="-clobber"
if test "$VISIT_VERSION" != current; then
  VISIT="$VISIT -v $VISIT_VERSION"
  XML2MAKEFILE_FLAGS="$XML2MAKEFILE_FLAGS -v $VISIT_VERSION"
fi
AC_SUBST(VISIT)

dnl ######################################################################
dnl
dnl Set private or public build
dnl
dnl ######################################################################

AC_ARG_ENABLE(public,
        AC_HELP_STRING([--enable-public],
        [Build the plugin in public use area]),
        XML2MAKEFILE_FLAGS="$XML2MAKEFILE_FLAGS -public"; ispublic=true, )
if test "$ispublic" = true; then
  case $host in
    nid*)
      AC_MSG_WARN(Make sure that your plugin is publicly readable.)
      ;;
  esac
  VS_TOPDIR=`dirname $VISIT_SYSDIR`
else
  VS_TOPDIR=$HOME/.visit
fi
AC_SUBST(VS_TOPDIR)

dnl ######################################################################
dnl
dnl Save variables and conclude
dnl
dnl ######################################################################

VISIT_DIR=`dirname $VISIT_BINDIR`
AC_SUBST(XML2MAKEFILE)
AC_SUBST(VISIT_DIR)
AC_SUBST(XML2MAKEFILE_FLAGS)
if test -n "$config_summary_file"; then
  echo "  VISIT:              $VISIT"              >> $config_summary_file
  echo "  XML2MAKEFILE:       $XML2MAKEFILE"       >> $config_summary_file
  echo "  XML2MAKEFILE_FLAGS: $XML2MAKEFILE_FLAGS" >> $config_summary_file
fi
