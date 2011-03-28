dnl###################################################################
dnl
dnl File:	tx_init.m4
dnl
dnl Purpose: Do standard initialization of autoconf projects including the
dnl          collection of data for keeping track of configuration metadata
dnl
dnl Version:	$Id: tx_init.m4 3654 2010-09-05 13:15:18Z cary $
dnl
dnl###################################################################

AC_PREFIX_DEFAULT(/usr/local/$PACKAGE_NAME)

dnl ######################################################################
dnl Set sourcedir and builddir for keeping track of the separation
dnl ######################################################################

abs_top_builddir=`pwd`
AC_SUBST(abs_top_builddir)
abs_top_srcdir0=`dirname $0`
abs_top_srcdir=`(cd $abs_top_srcdir0; pwd -P)`
AC_SUBST(abs_top_srcdir)
# Seems like abs_top_srcdir is being redefined by configure, so:
ABS_TOP_SRCDIR=$abs_top_srcdir
AC_SUBST(ABS_TOP_SRCDIR)
AC_DEFINE_UNQUOTED(ABS_TOP_SRCDIR, "$ABS_TOP_SRCDIR", [The top source directory])

dnl ######################################################################
dnl
dnl Get host
dnl
dnl ######################################################################

echo Getting canonical host
AC_CANONICAL_HOST
hosttype=$host
echo hosttype = $hosttype
AC_DEFINE_UNQUOTED(HOSTTYPE, "$host", [The host type])
echo HOSTTYPE = $HOSTTYPE
AC_MSG_CHECKING(hostname)
if hostnm=`/usr/ucb/hostname 2>/dev/null`; then
  :
elif hostnm=`hostname -f 2>/dev/null`; then
  :
else
  hostnm=`hostname`
fi
# At some computers, the above is not complete, so need to fix
case $hostnm in
  nid????? | hopper??)
    hostnm=$hostnm.nersc.gov
    ;;
esac
AC_MSG_RESULT($hostnm)
AC_DEFINE_UNQUOTED(HOSTNAME, "$hostnm", [The host name])
case $hostnm in
  *.intrepid.alcf.anl.gov)
    AC_DEFINE(IS_INTREPID, , [Whether this is intrepid.alcf.anl.gov])
    AM_CONDITIONAL(IS_INTREPID, true)
    ;;
  *)
    AM_CONDITIONAL(IS_INTREPID, false)
    ;;
esac

dnl ######################################################################
dnl
dnl Use wrapped automake
dnl
dnl ######################################################################

dnl # See whether we have automake
amver=`automake --version 2>/dev/null`
if test $? != 0; then
  echo automake not present in your path.
  echo Modifications to Makefile.am\'s will not propagate.
else
  AUTOMAKE=$abs_top_srcdir/config/automake.sh
fi
# For backware compatibility
EXEEXT=""
AC_SUBST(EXEEXT)

dnl ######################################################################
dnl
dnl Other nice initializations
dnl
dnl ######################################################################

dnl Prevent testing executables as fails on parallel machines
dnl JRC 20100107: Do we need this?  Probably needs to be added on
dnl a per-project basis
dnl AC_NO_EXECUTABLES

dnl ######################################################################
dnl
dnl Start configuration summary
dnl
dnl ######################################################################

builtin(include, config/tx_ac_summary.m4)
TX_AC_SUMMARY_INIT(config.summary, $0)

dnl ######################################################################
dnl
dnl Create configure string
dnl
dnl ######################################################################

if false; then
# Group args from quote to quote
composeargs() {
  unset composedargs
  while test -n "$1"; do
    newarg=`echo $1 | sed -e "s/^'//" -e "s/=/='/"`
    if echo $1 | grep "'\$" 1>/dev/null 2>&1; then
# Arg starts and ends with quote
      echo newarg = $newarg
      composedargs="$composedargs $newarg"
    else
      :
      # while test -n "$1"; do
      # shift
    fi
  done
}

composeargs $ac_configure_args

config_command="$0"
allargs=\"$ac_configure_args\"
echo allargs = $allargs
for i in $allargs; do
  arg="$i"
  echo arg = "$arg"
  if `echo $arg | grep ' ' 1>/dev/null 2>&1`; then
    if `echo $arg | grep = 1>/dev/null 2>&1`; then
      newarg=`echo $arg | sed -e "s/^'//" -e "s/=/='/"`
      echo newarg = $newarg
    fi
  else
    newarg=`echo $arg | sed -e "s/^'//" -e "s/'\$//"`
    echo newarg = $newarg
  fi
  config_command="$config_command $newarg"
done
echo config_command = $config_command
dnl AC_DEFINE_UNQUOTED(CONFIG_COMMAND, "$config_command", [The command used for configuration])
echo                          >>$config_summary_file
echo "Configuration command:" >>$config_summary_file
echo "$config_command"        >>$config_summary_file
fi

dnl ######################################################################
dnl
dnl Get search tools and version information
dnl
dnl ######################################################################

builtin(include, config/txsearch.m4)
#builtin(include, config/tx_svn_info.m4)
#TX_SVN_INFO([PROJECT_URL], [PROJECT_REV], .)
#tx_pkg_name=translit($PACKAGE,'a-z./-','A-Z____')
# eval "$tx_pkg_name"_URL=$PROJECT_URL
# eval "$tx_pkg_name"_REV=$PROJECT_REV
dnl AC_DEFINE_UNQUOTED(["$tx_pkg_name"_URL], "$PROJECT_URL", "SVN Project URL")
dnl AC_DEFINE_UNQUOTED(["$tx_pkg_name"_URL], "$PROJECT_REV", "SVN Project Revision")
#TX_SVN_INFO([CONFIG_URL], [CONFIG_REV], config)

dnl ######################################################################
dnl
dnl Supra search path
dnl
dnl ######################################################################

AC_MSG_CHECKING(SUPRA_SEARCH_PATH)
DEFAULT_SUPRA_SEARCH_PATH=$HOME/software:$HOME:/internal:/contrib:/usr/local:/opt/usr
if test -n "$UNIXFLAVOR"; then
  DEFAULT_SUPRA_SEARCH_PATH=$HOME/$UNIXFLAVOR:$DEFAULT_SUPRA_SEARCH_PATH
fi
AC_ARG_WITH(supra-search-path,
  AC_HELP_STRING([--with-supra-search-path=<supra-search-path>],
  [to look for installations under <supra-search-path>]),
  SUPRA_SEARCH_PATH="$withval")
if test -z "$SUPRA_SEARCH_PATH"; then
  SUPRA_SEARCH_PATH=$DEFAULT_SUPRA_SEARCH_PATH
fi
AC_ARG_WITH(extra-supra-search-path,
  AC_HELP_STRING([--with-extra-supra-search-path=<extra-supra-search-path>],
  [add extra path to default search path <extra-supra-search-path>]),
  EXTRA_SUPRA_SEARCH_PATH="$withval"
  )
if test -n "$EXTRA_SUPRA_SEARCH_PATH"; then
  SUPRA_SEARCH_PATH=$EXTRA_SUPRA_SEARCH_PATH:$SUPRA_SEARCH_PATH
fi
AC_MSG_RESULT($SUPRA_SEARCH_PATH)

echo >>$config_summary_file
echo Search for installations under the directories: >>$config_summary_file
TX_PRINT_VAR(SUPRA_SEARCH_PATH)


