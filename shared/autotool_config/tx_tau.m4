dnl ######################################################################
dnl
dnl Find tau wrappers
dnl
dnl ######################################################################

dnl Includes functions from txsearch.m4
builtin(include, config/txsearch.m4)

TAU_SP=$SUPRA_SEARCH_PATH
unset TAU_PATH
for i in `echo $TAU_SP | tr ':' ' '`; do
  TAU_PATH="$TAU_PATH:$i/tau"
done

unset TAU_BINDIR
AC_ARG_WITH(tau, AC_HELP_STRING([--with-tau],
	[Enable TAU performance measurement]),
	[ENABLE_TAU=yes], )

# If enabled, try to locate
if test "$ENABLE_TAU" = yes; then
  myproc=`uname -p`

  unset TAU_BINPATH
  for i in `echo $TAU_PATH | tr ':' ' '`; do
    TAU_BINPATH="$TAU_BINPATH:$i/bin"
  done

  # find tau_cxx.sh, the C++ compiler wrapper
  AC_PATH_PROG(TAU_CXX_WRAPPER, tau_cxx.sh, , $TAU_BINPATH)
  if test -z "$TAU_CXX_WRAPPER"; then
    AC_MSG_ERROR(tau_cxx.sh not found in $TAU_BINPATH.)
  fi

  # find tau_cc.sh, the C compiler wrapper
  AC_PATH_PROG(TAU_CC_WRAPPER, tau_cc.sh, , $TAU_BINPATH)
  if test -z "$TAU_CC_WRAPPER"; then
    AC_MSG_ERROR(tau_cc.sh not found in $TAU_BINPATH.)
  fi

  # find tau_f90.sh, the Fortran compiler wrapper
  AC_PATH_PROG(TAU_F90_WRAPPER, tau_f90.sh, , $TAU_BINPATH)
  if test -z "$TAU_F90_WRAPPER"; then
    AC_MSG_ERROR(tau_f90.sh not found in $TAU_BINPATH.)
  fi

  # set the TAU stub makefile
  TAU_DIR=${TAU_CXX_WRAPPER%/bin/tau_cxx.sh}
  TAU_MAKEFILE=${TAU_DIR}/include/Makefile
  if test ! -f $TAU_MAKEFILE; then
    AC_MSG_ERROR($TAU_MAKEFILE does not exist.)
  fi

  # set TAU instrumentation options
  tauopts="$tauopts -optVerbose -optPdtGnuFortranParser -optPreProcess -optNoCompInst"

  # replace CXX
  EXPCXX=`$TAU_CXX_WRAPPER -tau_makefile=$TAU_MAKEFILE -tau_options=\"$tauopts\" -show | sed 's/ .*$//'`
  AC_MSG_WARN(Tau expects you to be using $EXPCXX)
  CXX="$TAU_CXX_WRAPPER -tau_makefile=$TAU_MAKEFILE $tauopts $EXTRA_TAUOPTS $PROJECT_TAUOPTS"
  EXPCC=`$TAU_CC_WRAPPER -tau_makefile=$TAU_MAKEFILE -tau_options=\"$tauopts\" -show | sed 's/ .*$//'`
  AC_MSG_WARN(Tau expects you to be using $EXPCC)
  CC="$TAU_CC_WRAPPER -tau_makefile=$TAU_MAKEFILE $tauopts $EXTRA_TAUOPTS $PROJECT_TAUOPTS"
  AC_DEFINE(HAVE_TAU, , Define if have the Tau performance tools.)
  ABS_SERIALCC="$TAU_CC_WRAPPER -tau_makefile=$TAU_MAKEFILE $tauopts $EXTRA_TAUOPTS $PROJECT_TAUOPTS"
  ABS_SERIALCXX="$TAU_CXX_WRAPPER -tau_makefile=$TAU_MAKEFILE $tauopts $EXTRA_TAUOPTS $PROJECT_TAUOPTS"

  F77="$TAU_F90_WRAPPER -tau_makefile=$TAU_MAKEFILE $tauopts $EXTRA_TAUOPTS $PROJECT_TAUOPTS"
  FC="$TAU_F90_WRAPPER -tau_makefile=$TAU_MAKEFILE $tauopts $EXTRA_TAUOPTS $PROJECT_TAUOPTS"
  F90="$TAU_F90_WRAPPER -tau_makefile=$TAU_MAKEFILE $tauopts $EXTRA_TAUOPTS $PROJECT_TAUOPTS"
  FC_LD="$TAU_F90_WRAPPER -tau_makefile=$TAU_MAKEFILE $tauopts $EXTRA_TAUOPTS $PROJECT_TAUOPTS"

fi


AC_SUBST(TAU_CXX_WRAPPER)
AC_SUBST(TAU_MAKEFILE)
AM_CONDITIONAL(HAVE_TAU, test -n "$TAU_CXX_WRAPPER")

dnl ######################################################################
dnl
dnl Print out configuration information
dnl
dnl ######################################################################

if test -n "$config_summary_file"; then
   echo                                      >> $config_summary_file
   if test -n "$TAU_CXX_WRAPPER"; then
      echo "Using TAU with"                 >> $config_summary_file
      echo "  TAU_CXX_WRAPPER:    $TAU_CXX_WRAPPER"   >> $config_summary_file
      echo "  TAU_MAKEFILE:   $TAU_MAKEFILE"  >> $config_summary_file
      echo "  TAU_CXX:        $CXX"  >> $config_summary_file
   else
      echo "NOT using TAU"                >> $config_summary_file
   fi
fi

