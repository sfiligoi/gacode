dnl ######################################################################
dnl
dnl File:	    tx_babel.m4
dnl
dnl Purpose:	Determine where the babel files are.
dnl
dnl Version:	$Id: tx_babel.m4 3741 2010-11-19 20:30:41Z kruger $
dnl
dnl Copyright 2001-2010, Tech-X Corporation.  Redistribution allowed provided
dnl this copyright statement remains intact.
dnl
dnl ######################################################################

dnl default search paths for babel
# tx_babel_searchdir="$HOME"/software/babel:/contrib/babel:/usr/local/babel
BABEL_SHARED_SP=$SUPRA_SEARCH_PATH
unset tx_babel_searchdir
for i in `echo $BABEL_SHARED_SP | tr ':' ' '`; do
  tx_babel_searchdir="$tx_babel_searchdir:$i/babel-shared:$i/babel"
done

dnl AC_PATH_PROG(PYTHON, python)
dnl echo PYTHON = $PYTHON
PYDIR=python`python -c "import sys;print sys.version_info" | sed -e 's/^(//' -e 's/, /./' -e 's/,.*$//'`

dnl Locate all the babel related includes and libraries
TX_LOCATE_PKG(
  [babel-shared],
  [$tx_babel_searchdir],
  [sidl.h, sidl.hxx],
  [sidlstub_cxx, sidlstub_f90, sidl, chasmlite],
  [include/c:include/cxx:include/f77:include/f90:include/$PYDIR/babel],
  [])

dnl JRC: If Babel not found, nothing more to do except the AC_SUBST
if test -n "$BABEL_SHARED_LIBS"; then

dnl Look for sidl.mod in its various incarnations
  TX_PATH_FILES(BABEL_SHARED_SIDL_MOD, sidl.mod SIDL.mod, "", $BABEL_SHARED_INCDIR/f90:$BABEL_SHARED_INCDIR/../f90)
  if test -n "$BABEL_SHARED_SIDL_MOD"; then
    BABEL_SHARED_SIDL_MOD_DIR=`dirname $BABEL_SHARED_SIDL_MOD`
    if test $BABEL_SHARED_SIDL_MOD_DIR != $BABEL_SHARED_INCDIR; then
      BABEL_SHARED_INC="$BABEL_SHARED_INC -I$BABEL_SHARED_SIDL_MOD_DIR"
    fi
    AC_DEFINE(HAVE_SIDL_MOD, [], "Define if have sidl.mod file")
  fi
  if test `ls -d /usr/local/contrib/babel-shared/include/python*`; then
    BABELPYTHONDIR=`ls -d /usr/local/contrib/babel-shared/include/python*`
    BABEL_SHARED_INC="$BABEL_SHARED_INC -I$BABELPYTHONDIR/llnl_babel -I$BABELPYTHONDIR/llnl_babel_sidl_sidlx"
  fi
  if test -d "$PYTHON_LIBDIR/../site-packages/numpy/core/include/numpy"; then
    BABEL_SHARED_INC="$BABEL_SHARED_INC -I$PYTHON_LIBDIR/../site-packages/numpy/core/include/numpy"
  fi



dnl hacks to get this to work w/out having to add an extra
dnl configure option (--with-babel-bindir) if we want to use
dnl a 'special' babel.  You'd never want to use a babel binary
dnl that wasn't built the same way as lib. anyway.
dnl also changes so this works with distcomp2-SPM
  if test -n "$with_babel_dir"; then
    BABEL_SHARED_PATH=$with_babel_dir/bin
  else
    BABEL_SHARED_PATH=""
    for tx_itr in `echo "$tx_babel_searchdir" | tr ':' ' '`; do
      BABEL_SHARED_PATH=$BABEL_SHARED_PATH:${tx_itr}/bin
    done
  fi
  AC_PATH_PROGS(BABEL_SHARED_BIN, babel, , $BABEL_SHARED_PATH)

dnl duplicate for distcomp2 Makefile.am's -SPM
  BABEL_SHARED=$BABEL_SHARED_BIN
  BABEL_SHARED_BINDIR=`dirname "$BABEL_SHARED"`
  BABEL_SHARED_DIR=`dirname "$BABEL_SHARED_BINDIR"`
  BABEL_SHARED_INCDIR=$BABEL_SHARED_DIR/include
  BABEL_SHARED_LIBDIR=$BABEL_SHARED_DIR/lib
  BABEL_SHARED_LIBTOOL=$BABEL_SHARED_DIR/bin/babel-libtool
  BABEL_SHARED_CONFIG=$BABEL_SHARED_DIR/bin/babel-config
  BABEL_SHARED_CC=`$BABEL_SHARED_CONFIG --query-var=CC`
  BABEL_SHARED_CFLAGS=`$BABEL_SHARED_CONFIG --flags-c`
dnl JRC: error is below.  --version gives multiple lines.
dnl Must grep out the first.
dnl BABEL_SHARED_VERSION=`"$BABEL_SHARED" --version | cut -d ' ' -f3`
dnl Roopa: sed does not work in this case:
dnl        Babel version 1.4.0 (Revision: 6580 release-1-4-0-branch)
dnl BABEL_SHARED_VERSION=`"$BABEL_SHARED" --version | grep "^Babel version" | sed 's/^.* //'`
  BABEL_SHARED_VERSION=`"$BABEL_SHARED" --version | grep "^Babel version" | cut -d ' ' -f3`
  echo $BABEL_SHARED_VERSION
  BABEL_SHARED_LIB_SIDL_JAR=$BABEL_SHARED_LIBDIR/sidl-$BABEL_SHARED_VERSION.jar
  BABEL_SHARED_LIB_SIDLSTUB_JAR=$BABEL_SHARED_LIBDIR/sidlstub_$BABEL_SHARED_VERSION.jar
fi

AC_SUBST(BABEL_SHARED_BIN)
AC_SUBST(BABEL_SHARED)
AC_SUBST(BABEL_SHARED_BINDIR)
AC_SUBST(BABEL_SHARED_DIR)
AC_SUBST(BABEL_SHARED_INCDIR)
AC_SUBST(BABEL_SHARED_LIBDIR)
AC_SUBST(BABEL_SHARED_LIBTOOL)
AC_SUBST(BABEL_SHARED_CONFIG)
AC_SUBST(BABEL_SHARED_CC)
AC_SUBST(BABEL_SHARED_CFLAGS)
AC_SUBST(BABEL_SHARED_VERSION)
AC_SUBST(BABEL_SHARED_LIB_SIDL_JAR)
AC_SUBST(BABEL_SHARED_LIB_SIDLSTUB_JAR)
AM_CONDITIONAL(HAVE_BABEL_SHARED, test -n "$BABEL_SHARED")

dnl don't bother printing this if BABEL_SHARED was disabled
if test -n "$config_summary_file" -a x"$BABEL_SHARED_DOSEARCH" != xno; then
  printf "  Module files for babel:\n"          >> $config_summary_file
    TX_PRINT_VAR(BABEL_SHARED_SIDL_MOD)
    TX_PRINT_VAR(BABEL_SHARED_INC)
  printf "  Binaries sought: babel\n"          >> $config_summary_file
  if test -z "$BABEL_SHARED_BIN"; then
    AC_MSG_WARN([Unable to find shared babel binaries.  Use --with-babel-bindir to set the location.])
    TX_PRINT_VAR(BABEL_SHARED_BINDIR, -- failed --)
    TX_PRINT_VAR(BABEL_SHARED_BIN, -- failed --)
  else
    BABEL_SHARED_BINDIR=`dirname $BABEL_SHARED_BIN`
    TX_PRINT_VAR(BABEL_SHARED_BINDIR)
    TX_PRINT_VAR(BABEL_SHARED_BIN)
  fi

fi

dnl ##########################################################################
dnl
dnl Look for static babel
dnl
dnl ##########################################################################

BABEL_STATIC_SP=$SUPRA_SEARCH_PATH
unset tx_babel_static_searchdir
for i in `echo $BABEL_STATIC_SP | tr ':' ' '`; do
  tx_babel_static_searchdir="$tx_babel_static_searchdir:$i/babel-static"
done

dnl Locate all the babel related includes and libraries
TX_LOCATE_PKG(
  [babel-static],
  [$tx_babel_static_searchdir],
  [sidl.h, sidl.hxx],
  [sidlstub_cxx, sidlstub_f90, sidl, chasmlite],
  [include/c:include/cxx:include/f77:include/f90],
  [])

dnl we removed the sidl.mod test from above.  So to see if we add
dnl need the f90 directory add, we just test to see if it exists
dnl Note that the cxx directory can screw up the incdir
dnl if test -d "$BABEL_STATIC_INCDIR/f90" -o  -d "$BABEL_STATIC_INCDIR/../f90"; then
dnl   BABEL_STATIC_INC="$BABEL_STATIC_INC -I$BABEL_STATIC_INCDIR/f90"
dnl fi

dnl Find additional variables
if test -n "$BABEL_STATIC_LIBS"; then

dnl Look for sidl.mod in its various incarnations
  TX_PATH_FILES(BABEL_STATIC_SIDL_MOD, sidl.mod SIDL.mod, "", $BABEL_STATIC_INCDIR/f90:$BABEL_STATIC_INCDIR/../f90)
  if test -n "$BABEL_STATIC_SIDL_MOD"; then
    BABEL_STATIC_SIDL_MOD_DIR=`dirname $BABEL_STATIC_SIDL_MOD`
    if test $BABEL_STATIC_SIDL_MOD_DIR != $BABEL_STATIC_INCDIR; then
      BABEL_STATIC_INC="$BABEL_STATIC_INC -I$BABEL_STATIC_SIDL_MOD_DIR"
    fi
    AC_DEFINE(HAVE_STATIC_SIDL_MOD, [], "Define if have sidl.mod file for static babel")
  fi

dnl hacks to get this to work w/out having to add an extra
dnl configure option (--with-babel-bindir) if we want to use
dnl a 'special' babel.  You'd never want to use a babel binary
dnl that wasn't built the same way as lib. anyway.
dnl also changes so this works with distcomp2-SPM
  if test -n "$with_babel_static_dir"; then
    BABEL_STATIC_PATH=$with_babel_static_dir/bin
  else
    BABEL_STATIC_PATH=""
    for tx_itr in `echo "$tx_babel_static_searchdir" | tr ':' ' '`; do
      BABEL_STATIC_PATH=$BABEL_STATIC_PATH:${tx_itr}/bin
    done
  fi
  dnl echo BABEL_STATIC_LIBDIR = $BABEL_STATIC_LIBDIR
  BABEL_STATIC_BINDIR=`(cd $BABEL_STATIC_LIBDIR/../bin; pwd -P)`
  dnl echo BABEL_STATIC_BINDIR = $BABEL_STATIC_BINDIR
  BABEL_STATIC_PATH=$BABEL_STATIC_PATH:$BABEL_STATIC_BINDIR
  AC_PATH_PROGS(BABEL_STATIC_BIN, babel, , $BABEL_STATIC_PATH)
fi
if test -z "$BABEL_STATIC_BIN"; then
  BABEL_STATIC_BIN=`\ls $BABEL_STATIC_BINDIR/*-babel | sed 's/ .*$//'`
  if test -n "$BABEL_STATIC_BIN"; then
    AC_MSG_WARN(Is the Babel binary really $BABEL_STATIC_BIN?)
  fi
fi
# The lines below prevent debabelized build of FACETS, so I'm
# commenting them out. JAC 08/24/09
#if test -z "$BABEL_STATIC_BIN"; then
#  AC_MSG_ERROR(Could not find BABEL_STATIC_BIN)
#fi

dnl additional substitutions
AC_SUBST(BABEL_STATIC_BIN)
BABEL_STATIC_BINDIR=`dirname "$BABEL_STATIC_BIN"`
AC_SUBST(BABEL_STATIC_BINDIR)

printf "  Module files for babel:\n"          >> $config_summary_file
TX_PRINT_VAR(BABEL_STATIC_SIDL_MOD)
TX_PRINT_VAR(BABEL_STATIC_INC)

printf "  Binaries sought: babel\n"          >> $config_summary_file
if test -z "$BABEL_STATIC_BIN"; then
  AC_MSG_WARN([Unable to find babel binaries.  Use --with-babel-static-bindir to set the location.])
  TX_PRINT_VAR(BABEL_STATIC_BINDIR, -- failed --)
  TX_PRINT_VAR(BABEL_STATIC_BIN, -- failed --)
else
  TX_PRINT_VAR(BABEL_STATIC_BINDIR)
  TX_PRINT_VAR(BABEL_STATIC_BIN)
fi

dnl ##########################################################################
dnl
dnl CPP_BABEL needed in general
dnl
dnl ##########################################################################

if test -z "$CPP_BABEL"; then
  CPP_BABEL="gcc -E"
fi
AC_SUBST(CPP_BABEL)
TX_PRINT_VAR(CPP_BABEL)

dnl ##########################################################################
dnl
dnl generate babel implementation and client code during autotools build
dnl
dnl ##########################################################################

dnl
dnl top level checks.  don't try to do this if sidlDir not set
dnl enable-struct is temporary for distcomp2
dnl

if test -n "$BABEL_SHARED_BIN"; then
  if test ! -z ${sidlFile_0} ; then

    i=0
    while test $i -lt ${numGroup}; do
      sidlLocal=`eval 'echo $'"sidlFile_$i"`
      pathLocal=`eval 'echo $'"sidlPath_$i"`
      isCCA=`eval 'echo $'"useCCA_$i"`

      switch="-c -s"
      for sw in $switch; do

        if test "$sw" = "-c" ; then
          type="client"
          typeName=`eval 'echo $'"clientType_$i"`
        else
          type="impl"
          typeName=`eval 'echo $'"implType_$i"`
        fi

        dnl
        dnl parent directory of the sidl directory, where the generated
        dnl code will go
        dnl
        dirName=`dirname $sidlLocal`
        dirName=`dirname $dirName`
        dirName=${dirName}/`basename $sidlLocal .sidl`

        if test ! -z $pathLocal; then
          dirName=`eval 'echo $'"pathLocal"`
        fi
        dnl
        dnl run babel
        dnl
        if test "$isCCA" = "yes" ; then
          echo "${CCAFE_BINDIR}/babel ${sw}${typeName} -R${CCASPEC_BABEL_XML_REPOSITORY} -o ${abs_top_builddir}/$dirName/$type/$typeName ${abs_top_srcdir}/${sidlLocal}"
          ${CCAFE_BINDIR}/babel ${sw}${typeName} -R${CCASPEC_BABEL_XML_REPOSITORY} -o ${abs_top_builddir}/$dirName/$type/$typeName ${abs_top_srcdir}/${sidlLocal}
        else
          echo "${BABEL_SHARED_BINDIR}/babel ${sw}${typeName} -o ${abs_top_srcdir}/$dirName/$type/$typeName ${abs_top_srcdir}/${sidlLocal}"
          ${BABEL_SHARED_BINDIR}/babel ${sw}${typeName} -o ${abs_top_srcdir}/$dirName/$type/$typeName ${abs_top_srcdir}/${sidlLocal}
        fi
        cd ${abs_top_builddir}
      done
      dnl loop over sidl directories
      i=`expr $i + 1`
    done
    dnl loop over binding type

  fi
fi

dnl end top level check

