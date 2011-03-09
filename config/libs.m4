dnl ######################################################################
dnl
dnl File:	libs.m4
dnl
dnl Purpose:	Determine how to build libraries, esp. shared
dnl
dnl Version:	$Id: libs.m4 3780 2011-01-21 22:00:05Z srinath $
dnl
dnl Copyright 2003-2010, Tech-X Corporation.  This file may be freely
dnl redistributed and modified provided the above copyright remains.
dnl
dnl ######################################################################

#
# Determine whether to use rpath
#
AC_ARG_ENABLE(rpath, AC_HELP_STRING([--disable-rpath],
  [to not link with a runtime link path]),
  [USE_RPATH="${enableval}"], USE_RPATH=yes)

#
# Global find of some files
#

# Get any user specified variables
AC_ARG_WITH(z-dir, [  --with-z-dir=<zlib installation directory>],
        Z_DIR="$withval")
AC_ARG_WITH(z-libdir, [  --with-z-libdir=<z library directory>],
        Z_LIBDIR="$withval")
if test -n "$Z_LIBDIR"; then
  Z_LIBPATH=$Z_LIBDIR
elif test -n "$Z_DIR"; then
  Z_LIBPATH=$Z_DIR/lib
else
  Z_LIBPATH=/usr/lib:/usr/local/lib:/usr/common/usg/gnu/lib:/usr/lib64
fi

builtin(include, config/txsearch.m4)
TX_PATH_FILES(LIBZ, libz.a libz.so libz.dylib, "", $Z_LIBPATH)
if false; then
if test -z "$LIBZ"; then
  TX_PATH_FILES(LIBZ, libz.so, "", $Z_LIBPATH)
fi
if test -z "$LIBZ"; then
  TX_PATH_FILES(LIBZ, libz.dylib, "", $Z_LIBPATH)
fi
fi
if test -n "$LIBZ"; then
  Z_LIB=-lz
  Z_LIBDIR=`dirname $LIBZ`
  Z_LIBS="-L$Z_LIBDIR -lz"
fi
TX_CLEAN_LIBS([Z_LIBS])

# Libraries specific to MinGW, e.g. -lwinmm
MINGW_LIBS=""

# Initial blank value for SEARCH_ORDER
SEARCH_ORDER=""

# Needed by automake 1.5.  Should be done on
# individual basis
OBJEXT=o

case "$host" in
  xt3-*)
    NET_LIBS=""
    THREAD_LIB=""
    DB_LIB=""
    DYNLINK_LIB=""
    LINK_STATIC=""
    LINK_DYNAMIC=""
    RPATH_FLAG="-L"
    HAVE_BUNDLES=false
    LINKCMODULE="$CC"
    LINKCMODULE_FLAGS=""
    LINKCXXMODULE="$CXX"
    LINKCXXMODULE_FLAGS=""
    LDWITHMODULE="$CXX"
    LDWITHMODULE_FLAGS=""
    ;;

  alpha*-dec-osf*)
    NET_LIBS="-lbsd"
    THREAD_LIB="-lpthread"
    DB_LIB="-ldb"
    DYNLINK_LIB=""
    SO=".so"
    HAVE_BUNDLES=false

    case $SERIALCXX in
      cxx)
	PIC_FLAG=" "
        RPATH_FLAG="-rpath "
        LINKCXXMODULE="$CXX -shared -expect_unresolved \"*\""
        LDWITHMODULE="$CXX "
        LDWITHMODULE_FLAGS=" "
	;;
      g++)
	PIC_FLAG="-fPIC"
        RPATH_FLAG=-Wl,-rpath,
        LINKCXXMODULE="$CXX -shared -Wl,-expect_unresolved,\"*\""
        LINKCXXMODULE_FLAGS="-shared -Wl,-expect_unresolved,\"*\""
        LDWITHMODULE="$CXX"
        LDWITHMODULE_FLAGS=" "
	;;
      KCC)
	PIC_FLAG=" "
        RPATH_FLAG="-rpath "
        LINKCXXMODULE="$CXX -shared "
        LDWITHMODULE="$CXX "
        LDWITHMODULE_FLAGS=" "
	;;
    esac

    case $SERIALCC in
      cc | */cc)
        LINKCMODULE="$CC -shared -expect_unresolved \"*\""
	;;
      gcc | */gcc)
        LINKCMODULE="$CC -shared -Wl,-expect_unresolved,\"*\""
        LINKCMODULE_FLAGS="-shared -Wl,-expect_unresolved,\"*\""
	;;
    esac
    ;;

  *-darwin*)
    NET_LIBS=""
    THREAD_LIB="-lpthread"
    DB_LIB="-ldb"
    DYNLINK_LIB="-ldl"
    SHAREDLIBFLAGS="-dynamiclib"
    SO=".dylib"
    LINK_STATIC="-Wl,-static"
    LINK_DYNAMIC="-Wl,-dynamic"
    LDUNDEFSYM="-undefined dynamic_lookup"
    # RPATH_FLAG="-L"
    RPATH_FLAG="-Wl,-rpath,"	# As of 10.5
    HAVE_BUNDLES=true
    SEARCH_ORDER="-search_paths_first"
    case $SERIALCXX in
      *g++ | *c++ | *g++3)
	PIC_FLAG="-fpic"
        LINKCMODULE="$CC -shared "
        LINKCMODULE_FLAGS="-shared "
        LINKCXXMODULE="$CXX -shared "
        LINKCXXMODULE_FLAGS="-shared "
# The following flags are unknown as of 12 July 02
        # LDWITHMODULE="$CXX -Wl,-export-dynamic"
        # LDWITHMODULE_FLAGS="-Wl,-export-dynamic"
        LDWITHMODULE="$CXX "
        LDWITHMODULE_FLAGS=""
	;;
    esac
    ;;

  *-hp-hpux10*)
    NET_LIBS="-lnsl -ldld"
    THREAD_LIB="-lpthread"
    DB_LIB="-ldbm"
    DYNLINK_LIB="-lnsl -ldld"
    SO=".sl"
    HAVE_BUNDLES=false

    case $SERIALCXX in
      aCC)
        ver=`aCC -V 2>&1 | sed "s/^.*A\.//"`
        minver=`echo $ver | sed "s/^.*\.//"`
        majver=`echo $ver | sed "s/\..*$//"`
        if test $majver -lt 3 -o $majver -eq 3 -a $minver -lt 13; then
          AC_MSG_WARN(aCC not tested for versions less than A.03.13)
        fi
	;;
      g++)
	PIC_FLAG="-fpic"
	RPATH_FLAG=-Wl,+b,
        gccspecdir=`$CC -v 2>&1 | grep specs | sed -e "s/^.* //"`
        gcclibdir=`echo $gccspecdir | sed -e "s/\/specs//"`
	CXXMODULELIBS="-L$gcclibdir -lstdc++"
        LINKCXXMODULE="ld -b"
        LDWITHMODULE="$CXX -Wl,-E -Wl,+s"
        LDWITHMODULE_FLAGS="-Wl,-E -Wl,+s"
	;;
    esac
    ;;

  *-hp-hpux11*)
    NET_LIBS="-lnet"
    THREAD_LIB=" "
    # No working thread lib on hp?
    DB_LIB="-ldbm"
    DYNLINK_LIB="-lnsl -ldld"
    SO=".sl"
    HAVE_BUNDLES=false

    case $SERIALCXX in
      aCC)
        ver=`aCC -V 2>&1 | sed "s/^.*A\.//"`
        minver=`echo $ver | sed "s/^.*\.//"`
        majver=`echo $ver | sed "s/\..*$//"`
        if test $majver -lt 3 -o $majver -eq 3 -a $minver -lt 13; then
          AC_MSG_WARN(aCC not tested for versions less than A.03.13)
        fi
	PIC_FLAG="+z"
	RPATH_FLAG="+b "
        LINKCXXMODULE="$CXX -b"
	CXXMODULELIBS="-lm"
        LDWITHMODULE="$CXX -E"
        LDWITHMODULE_FLAGS="-E"
	;;
      g++)
	AC_MSG_ERROR($SERIALCXX not supported on $host)
	;;
    esac

    case $SERIALCC in
      acc | */acc)
        LINKCMODULE="$CC -b"
	CMODULELIBS="-lm"
	;;
      gcc | */gcc)
        # LINKCMODULE="$CC -Wl,-b"
	# Above fails trying to link in position dependent code in /usr/ccs/lib/crt0.o
        gccspecdir=`$CC -v 2>&1 | grep specs | sed -e "s/^.* //"`
        gcclibdir=`echo $gccspecdir | sed -e "s/\/specs//"`
	CMODULELIBS="-L$gcclibdir -lgcc -lm"
        LINKCMODULE="ld -b"
	;;
    esac
    ;;

  *-ibm-aix*)
    RPATH_FLAG=-L
    PIC_FLAG=""
    NET_LIBS=" "
    THREAD_LIB="-lpthreads"
    DB_LIB="-ldb"
    DYNLINK_LIB="-lld -ldl"
    SO=".so"
    HAVE_BUNDLES=false
    # echo SERIALCXX = $SERIALCXX
    case $SERIALCXX in

      *g++)
	LDFLAGS="$LDFLAGS -Wl,-bbigtoc"
        LINKCMODULE='$(PYTHON_LIBDIR)/ld_so_aix $(CC) -bI:$(PYTHON_LIBDIR)/python.exp'
        LINKCXXMODULE='$(PYTHON_LIBDIR)/ld_so_aix $(CXX) -bI:$(PYTHON_LIBDIR)/python.exp'
        LDWITHMODULE='$(PYTHON_LIBDIR)/makexp_aix python.exp "" $(PYTHON_LIBDIR)/lib$(PYTHON_LIB).a; $(CXX) -Wl,-bE:python.exp'
	;;

      *xlC)
	AC_MSG_WARN(You must add the -C option to nm in $PYTHON_LIBDIR/makexp_aix)
        LINKCMODULE='$(PYTHON_LIBDIR)/ld_so_aix $(CC) -bI:$(PYTHON_LIBDIR)/python.exp'
	# Note that we use C compiler, not C++ compiler, below
        LINKCXXMODULE='$(PYTHON_LIBDIR)/ld_so_aix $(CXX) -bI:$(PYTHON_LIBDIR)/python.exp'
        LDWITHMODULE='$(PYTHON_LIBDIR)/makexp_aix python.exp "" $(PYTHON_LIBDIR)/lib$(PYTHON_LIB).a; $(CXX) -bE:python.exp'
	;;

      *KCC)
	;;

    esac
    ;;

  *-cygwin*)
    RPATH_FLAG="-Wl,-rpath,"
    DYNLINK_LIB="-ldl"
    ;;

  *-pc-mingw*)
    RPATH_FLAG="-Wl,-rpath,"
    DYNLINK_LIB=""
    MINGW_LIBS="-lwinmm -lws2_32"
    ;;

  *-linux* | *-freebsd*)
    linuxver=`uname -r | sed 's/-.*$//' | sed 's/\.[0-9]*$//'`
    linuxminver=`echo $linuxver | sed 's/^.*\.//'`
    linuxmajver=`echo $linuxver | sed 's/\..*$//'`
    NET_LIBS=""
    if test "$linuxmajver" -eq 2 -a "$linuxminver" -le 2; then
      NET_LIBS="-lbsd"
    fi
    THREAD_LIB="-lpthread"
    DB_LIB="-ldb"
    DYNLINK_LIB="-ldl"
    SHAREDLIBFLAGS="-shared"
    SO=".so"
    LINK_STATIC="-Wl,-Bstatic"
    LINK_DYNAMIC="-Wl,-Bdynamic"
    LDUNDEFSYM="--no-undefined"
    RPATH_FLAG="-Wl,-rpath,"
    HAVE_BUNDLES=false
    case $SERIALCXX in
      *g++ | *c++)
	PIC_FLAG="-fpic"
	RPATH_FLAG="-Wl,-rpath,"
        LINKCMODULE="$CC -shared "
        LINKCMODULE_FLAGS="-shared "
        LINKCXXMODULE="$CXX -shared "
        LINKCXXMODULE_FLAGS="-shared "
        LDWITHMODULE="$CXX -Wl,-export-dynamic"
        LDWITHMODULE_FLAGS="-Wl,-export-dynamic"
	;;
    esac
    ;;

  *-sgi-irix6*)
    NET_LIBS="-lbsd"
    THREAD_LIB="-lpthread"
    DB_LIB="-ldb"
    DYNLINK_LIB=""
    SO=".so"
    HAVE_BUNDLES=false

    case $SERIALCXX in
      CC)
	PIC_FLAG=" "
	RPATH_FLAG="-rpath "
        LINKCMODULE="$CC -shared"
        LINKCXXMODULE="$CXX -shared"
        LDWITHMODULE="$CXX"
        LDWITHMODULE_FLAGS=" "
	;;
      g++)
	PIC_FLAG="-fpic"
	RPATH_FLAG=-Wl,-rpath,
        LINKCMODULE="$CC -shared"
        LINKCXXMODULE="$CXX -shared"
        LDWITHMODULE="$CXX "
        LDWITHMODULE_FLAGS=" "
	;;
    esac
    ;;

  *-*-solaris2.5* | *-*-solaris2.6*)
    NET_LIBS="-lsocket -lnsl"
    THREAD_LIB="-lpthread"
    DB_LIB="-ldb"
    DYNLINK_LIB="-ldl"
    SO=".so"
    HAVE_BUNDLES=false

    case $SERIALCXX in
      CC)
	PIC_FLAG="-PIC"
	RPATH_FLAG=-Wl,-R,
        LINKCXXMODULE="$CXX -G"
        LDWITHMODULE="$CXX"
        LDWITHMODULE_FLAGS=" "
	;;
      g++)
	PIC_FLAG="-fpic"
	RPATH_FLAG=-Wl,-R,
	LDFLAGS="$LDFLAGS -Wl,-z,muldefs "
        LINKCXXMODULE="$CXX -G -nostartfiles -Wl,-z,muldefs"
        LDWITHMODULE="$CXX -Wl,-z,muldefs"
        LDWITHMODULE_FLAGS="-Wl,-z,muldefs"
	;;
    esac

    case $SERIALCC in
      cc | */cc)
        LINKCMODULE="$CC -G"
	;;
      gcc | */gcc)
        LINKCMODULE="$CC -G -nostartfiles -Wl,-z,muldefs"
	;;
    esac
    ;;

  *-*-solaris2.7* |  *-*-solaris2.8*)
    NET_LIBS="-lsocket -lnsl"
    THREAD_LIB="-lpthread"
    DB_LIB="-ldb"
    DYNLINK_LIB="-ldl"
    SO=".so"
    HAVE_BUNDLES=false

    case $SERIALCXX in
      CC)
	PIC_FLAG="-PIC"
	RPATH_FLAG=-Wl,-R,
        LINKCXXMODULE="$CXX -G"
        LDWITHMODULE="$CXX"
        LDWITHMODULE_FLAGS=" "
	;;
      g++)
	PIC_FLAG="-fpic"
	RPATH_FLAG=-Wl,-R,
        LINKCXXMODULE="$CXX -G -nostartfiles -Wl,-z,muldefs"
        # LDWITHMODULE="$CXX -Wl,-z,muldefs -Wl,-export-dynamic"
        # LDWITHMODULE_FLAGS="-Wl,-z,muldefs -Wl,-export-dynamic"
        LDWITHMODULE="$CXX -Wl,-z,muldefs"
        LDWITHMODULE_FLAGS="-Wl,-z,muldefs"
	;;
    esac

    case $SERIALCC in
      cc | */cc)
        LINKCMODULE="$CC -G"
	;;
      gcc | */gcc)
        LINKCMODULE="$CC -G -nostartfiles -Wl,-z,muldefs"
	;;
    esac
    ;;

  *-cray-unicos*)
    NET_LIBS=""
    THREAD_LIB="-lpthread"
    DB_LIB="-ldb"
    DYNLINK_LIB=""
    SO=".so"
    PIC_FLAG="-fpic"
    LINKCMODULE="$CC -shared "
    LINKCMODULE_FLAGS="-shared "
    LINKCXXMODULE="$CXX -shared "
    LINKCXXMODULE_FLAGS="-shared "
    LDWITHMODULE="$CXX -Wl,-export-dynamic"
    LDWITHMODULE_FLAGS="-Wl,-export-dynamic"
    HAVE_BUNDLES=false
    case $SERIALCXX in
      KCC)
        RPATH_FLAG="-L"
	;;
      CC)
        RPATH_FLAG="-L"
	;;
    esac
    ;;

  *)
    AC_MSG_WARN(Libraries unknown for host $host.  Please notify the developers.)
    ;;

esac

# For g++, add in location of C++ libraries
case $SERIALCXX in
  *g++)
    gcclibplatform=`$CXX -dumpmachine`          #osx 10.4 + gcc4.0.1 changed libstdc++.a to libstdc++-static.a
    gcclibversion=`$CXX -dumpversion`
    mysystype="${gcclibplatform}-${gcclibversion}"
    case ${mysystype} in
      *-apple-darwin8-4.0.*)
        gcclibextention='-static'
        echo 'mac-gcc 4.0.1: we are looking for libstdc++-static.a'
        ;;
      *)
        gcclibextention=''
        ;;
    esac
    gcclibfilename=`$CXX -print-file-name=libstdc++${gcclibextention}.a`
    # echo gcclibfilename = $gcclibfilename
    gcclibsdir=`dirname $gcclibfilename`
    gcclibsdir=`(cd $gcclibsdir; pwd -P)`
    # echo gcclibsdir = $gcclibsdir
# Need the -L flag or the system library might be picked up
# instead of a specially built library.
    COMPILER_LIBDIR=${gcclibsdir}
    COMPILER_LIBFLAG="${RPATH_FLAG}${gcclibsdir} -L${gcclibsdir}"
    COMPILER_LTLIBFLAG="-rpath ${gcclibsdir}"
# For backward compatibility
    CXX_LIBDIR=$COMPILER_LIBDIR
    CXX_LIBFLAG="$COMPILER_LIBFLAG"
    CXX_LTLIBFLAG="$COMPILER_LTLIBFLAG"
    ;;

esac


# In case wrapper compiler uses shared libraries without setting rpath:
if test -n "$MPI_LIBDIR"; then
  MPI_RUNLIBFLAG=${RPATH_FLAG}$MPI_LIBDIR
  LDFLAGS="$LDFLAGS $MPI_RUNLIBFLAG"
fi
# JRC: For backward compatibility
ZLIB_DIR=$Z_LIBDIR

# Put into cache
AC_SUBST(CMODULELIBS)
AC_SUBST(COMPILER_LIBFLAG)
AC_SUBST(COMPILER_LTLIBFLAG)
AC_SUBST(CXXMODULELIBS)
AC_SUBST(CXX_LIBDIR)
AC_SUBST(CXX_LIBFLAG)
AC_SUBST(CXX_LTLIBFLAG)
AC_SUBST(DB_LIB)
AC_SUBST(DYNLINK_LIB)
AC_SUBST(LDUNDEFSYM)
AC_SUBST(HAVE_BUNDLES)
AC_SUBST(LDWITHMODULE)
AC_SUBST(LDWITHMODULE_FLAGS)
AC_SUBST(LINKCMODULE)
AC_SUBST(LINKCMODULE_FLAGS)
AC_SUBST(LINKCXXMODULE)
AC_SUBST(LINKCXXMODULE_FLAGS)
AC_SUBST(LINK_DYNAMIC)
AC_SUBST(LINK_STATIC)
AC_SUBST(MINGW_LIBS)
AC_SUBST(MPI_RUNLIBFLAG)
AC_SUBST(NET_LIBS)
AC_SUBST(OBJEXT)
AC_SUBST(PIC_FLAG)
AC_SUBST(RPATH_FLAG)
AC_SUBST(SEARCH_ORDER)
AC_SUBST(SO)
AC_SUBST(THREAD_LIB)
AC_SUBST(ZLIB_DIR)
AC_SUBST(Z_LIBDIR)
AC_SUBST(Z_LIBS)

# Print results
echo >>$config_summary_file
echo "Various link flags" >>$config_summary_file
TX_PRINT_VAR(CMODULELIBS)
TX_PRINT_VAR(COMPILER_LIBFLAG)
TX_PRINT_VAR(COMPILER_LTLIBFLAG)
TX_PRINT_VAR(CXXMODULELIBS)
TX_PRINT_VAR(CXX_LIBDIR)
TX_PRINT_VAR(CXX_LIBFLAG)
TX_PRINT_VAR(CXX_LTLIBFLAG)
TX_PRINT_VAR(DB_LIB)
TX_PRINT_VAR(DYNLINK_LIB)
TX_PRINT_VAR(LDUNDEFSYM)
TX_PRINT_VAR(HAVE_BUNDLES)
TX_PRINT_VAR(LDWITHMODULE)
TX_PRINT_VAR(LDWITHMODULE_FLAGS)
TX_PRINT_VAR(LINKCMODULE)
TX_PRINT_VAR(LINKCMODULE_FLAGS)
TX_PRINT_VAR(LINKCXXMODULE)
TX_PRINT_VAR(LINKCXXMODULE_FLAGS)
TX_PRINT_VAR(LINK_DYNAMIC)
TX_PRINT_VAR(LINK_STATIC)
TX_PRINT_VAR(MINGW_LIBS)
TX_PRINT_VAR(MPI_RUNLIBFLAG)
TX_PRINT_VAR(NET_LIBS)
TX_PRINT_VAR(OBJEXT)
TX_PRINT_VAR(PIC_FLAG)
TX_PRINT_VAR(RPATH_FLAG)
TX_PRINT_VAR(SEARCH_ORDER)
TX_PRINT_VAR(SO)
TX_PRINT_VAR(THREAD_LIB)
TX_PRINT_VAR(ZLIB_DIR)
TX_PRINT_VAR(Z_LIBDIR)
TX_PRINT_VAR(Z_LIBS)
TX_PRINT_VAR(Z_RPLIBS)
TX_PRINT_VAR(Z_LTLIBS)
TX_PRINT_VAR(Z_ALIBS)

dnl ######################################################################
dnl
dnl Get MPI libraries
dnl
dnl ######################################################################

REMOVE_CXX_MPI_LIBS=false
REMOVE_FC_MPI_LIBS=false
if test "$parallel" = yes; then

  unset CXX_MPI_LIBS
  case $CXXBASE in
    mp*)
      cands=`$CXX -show 2>/dev/null`
      if test $? = 0; then
        echo MPI candidates are $cands
        for i in $cands; do
          # echo  Examining $i
          case $i in
            -L* | -l* | /*.a)
              CXX_MPI_LIBS="$CXX_MPI_LIBS $i"
              ;;
          esac
        done
      fi
      REMOVE_CXX_MPI_LIBS=true
      ;;
  esac
  TX_CLEAN_LIBS([CXX_MPI_LIBS])
# Extra substitution does not hurt
  if test -f $config_summary_file; then
    echo >>$config_summary_file
    echo CXX_MPI_LIBS and derivatives >>$config_summary_file
    TX_PRINT_VAR(CXX_MPI_LIBS)
    TX_PRINT_VAR(CXX_MPI_RPLIBS)
    TX_PRINT_VAR(CXX_MPI_LTLIBS)
    TX_PRINT_VAR(CXX_MPI_ALIBS)
  fi

  unset FC_MPI_LIBS
  case $FCBASE in
    mp*)
      cands=`$FC -show 2>/dev/null`
      if test $? = 0; then
        echo MPI candidates are $cands
        for i in $cands; do
          # echo  Examining $i
          case $i in
            -L* | -l* | /*.a)
              FC_MPI_LIBS="$FC_MPI_LIBS $i"
              ;;
          esac
        done
      fi
      REMOVE_FC_MPI_LIBS=true
      ;;
  esac
  TX_CLEAN_LIBS([FC_MPI_LIBS])
# Extra substitution does not hurt
  if test -f $config_summary_file; then
    echo >>$config_summary_file
    echo FC_MPI_LIBS and derivatives >>$config_summary_file
    TX_PRINT_VAR(FC_MPI_LIBS)
    TX_PRINT_VAR(FC_MPI_RPLIBS)
    TX_PRINT_VAR(FC_MPI_LTLIBS)
    TX_PRINT_VAR(FC_MPI_ALIBS)
  fi

fi

AC_SUBST([CXX_MPI_LIBS])
AC_SUBST([CXX_MPI_RPLIBS])
AC_SUBST([CXX_MPI_LTLIBS])
AC_SUBST([CXX_MPI_ALIBS])

AC_SUBST([FC_MPI_LIBS])
AC_SUBST([FC_MPI_RPLIBS])
AC_SUBST([FC_MPI_LTLIBS])
AC_SUBST([FC_MPI_ALIBS])

