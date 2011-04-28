dnl ######################################################################
dnl
dnl
dnl File:	corba.m4
dnl
dnl Purpose:	Find CORBA stuff for Makefiles
dnl
dnl Version:	$Id: corba.m4 3366 2010-01-15 18:43:10Z dws $
dnl
dnl Copyright 2001-2010, Tech-X Corporation.  This file may be freely
dnl distributed provided copyright statement remains in place.
dnl
dnl ######################################################################

dnl ######################################################################
dnl
dnl Get the IDL to C++ translator
dnl
dnl ######################################################################

builtin(include, config/txsearch.m4)

AC_MSG_WARN(If this hangs - you may have RSI's idl (Interactive Data Language) in your path ahead of the Interface Definition Language translator for CORBA.  You will need to change your path.)

##########
#
# CORBA_IDL_PROG is the name of the idl translator.  It can be
# set to tao_idl to require use of TAO
#
##########

CORBA_TAO_IDL_PROG=tao_idl
if test -z "$CORBA_IDL_PROG"; then
  CORBA_IDL_PROG=idl
fi

##########
#
# Allow user to set absolute path to the idl translator.
# Otherwise look in standard places
#
##########

AC_ARG_WITH(corba-idl,
[  --with-corba-idl=<corba-idl>  set location of the CORBA idl translator],
CORBA_IDL="$withval")
if test -n "$CORBA_IDL"; then
  IDL=$CORBA_IDL
else
  if test -n "$USEOB3"; then
    IDL_PATH=$HOME/$host/OB3-${COMPDIR}/bin:$HOME/$host/OB3/bin:/usr/local/OB3-${COMPDIR}/bin:/usr/local/OB3/bin:/loc/OB3-${COMPDIR}/bin:/local/OB3-${COMPDIR}/bin:/loc/OB3/bin:/mfe/local/OB3-${COMPDIR}/bin:/mfe/local/OB3/bin:$PATH
  else
    IDL_PATH=$HOME/$host/OB4-${COMPDIR}/bin:$HOME/$host/OB4/bin:/usr/local/OB4-${COMPDIR}/bin:/usr/local/OB4/bin:/usr/local/ACE+TAO/bin:/local/OB4-${COMPDIR}/bin:/loc/OB4-${COMPDIR}/bin:/loc/OB4/bin:/mfe/local/OB4-${COMPDIR}/bin:/mfe/local/OB4/bin:$PATH
  fi
  AC_PATH_PROGS(IDL, $CORBA_IDL_PROG, "", $IDL_PATH)
  if test -z "$IDL"; then
    AC_PATH_PROGS(IDL, $CORBA_TAO_IDL_PROG, "", $IDL_PATH)
  fi
fi

IDL_DIR=`dirname $IDL`

# Try to get version
# Test whether omniORB, TAO or MICO is being used
isomni=`echo $IDL | grep omni`
istao=`echo $IDL | grep tao`
ismico=`echo $IDL | grep mico`

# Distinguish among various ORBs
if test -n "$ismico"; then
  echo ... Using the MICO IDL to C++ Translator
  CXXORB=MICO
  TESTINC=CORBA.h
else
  if test -n "$isomni"; then
    echo ... Using the omniORB IDL to C++ Translator
    CXXORB=omniORB
    TESTINC=CORBA.h
  else
    if test -n "$istao"; then
      echo ... Using the TAO IDL to C++ Translator
      CXXORB=tao
      TESTINC=corba.h

      # Truncate IDL dir again to get base dir
      BASE_DIR=`dirname $IDL_DIR`

      if test "x$ACE_ROOT" = "x"; then
        if test -e "$BASE_DIR/include/ace/ACE.h" -o -e "$BASE_DIR/ace/ACE.h"; then
          ACE_ROOT=$BASE_DIR
        else
	  AC_MSG_ERROR([Could not find valid ACE_ROOT for tao_idl $IDL with base directory $BASE_DIR])
        fi
      else
        if test "$ACE_ROOT" != "$BASE_DIR"; then
	   AC_MSG_ERROR([ACE_ROOT($ACE_ROOT) does not match $BASE_DIR for $IDL])
	fi
      fi

      if test "x$TAO_ROOT" = "x"; then
        AC_MSG_NOTICE([Environment variable TAO_ROOT is not defined.])
	if test -d $ACE_ROOT/TAO; then
           TAO_ROOT=$ACE_ROOT/TAO
	   AC_MSG_NOTICE([... Standard build, setting TAO_ROOT=$TAO_ROOT])
	else
	   TAO_ROOT=$ACE_ROOT
	   AC_MSG_NOTICE([... Could be autotools build, setting TAO_ROOT=$TAO_ROOT])
	fi
      elif test "$ACE_ROOT/TAO" != "$TAO_ROOT"; then
      	  # Maybe autotool install
	  if test "$ACE_ROOT" != "$TAO_ROOT"; then
	     AC_MSG_ERROR([Environment variable TAO_ROOT($TAO_ROOT) not valid for ACE_ROOT($ACE_ROOT)])
	  fi
      fi

      echo TAO_ROOT=$TAO_ROOT
      echo ACE_ROOT=$ACE_ROOT

      # Determine TAO LD-PATH
      TX_PATH_FILES(ACE_LIB_PATH, libTAO.a libTAO.so libTAO.dylib, "", "$ACE_ROOT/lib:TAO_ROOT/lib")
      if test -z "$ACE_LIB_PATH"; then
        AC_MSG_ERROR([Cannot locate TAO libraries])
      fi
      ACE_LIB_PATH=`dirname $ACE_LIB_PATH`

      # Preset platform-specific path

      case "$host" in
        *-*-darwin*)
          IDL="ACE_ROOT=$ACE_ROOT TAO_ROOT=$TAO_ROOT DYLD_LIBRARY_PATH=$IDL_DIR:$ACE_LIB_PATH:$DYLD_LIBRARY_PATH $IDL"
          ;;
        i686-*-linux* | i386-*-linux*)
          IDL="ACE_ROOT=$ACE_ROOT TAO_ROOT=$TAO_ROOT LD_LIBRARY_PATH=$ACE_LIB_PATH:$LD_LIBRARY_PATH $IDL "
          case "$SERIALCXX" in
            *++ | *g++3)
              CORBA_CXXFLAGS="-D_POSIX_THREAD -D_POSIX_THREAD_SAFE_FUNCTIONS -D_REENTRANT -DACE_HAS_AIO_CALLS"
              ;;
           esac
           ;;
        *-sgi-irix6*)
          # IDL="LD_LIBRARY_PATH=$ACE_LIB_PATH/ace:$LD_LIBRARY_PATH $IDL -Ge 1 -Sc"
          CORBA_CXXFLAGS="-DACE_HAS_EXCEPTIONS -diag_suppress 3284 -ptused -prelink"
          ;;
        *)
          ;;
      esac

    else
      # Check if we are ORBacus
      IDLVER=`$IDL --version 2>&1`

      isorbacus=`echo $IDLVER | grep ORBacus`
      if test -n "$isorbacus"; then
        echo ... Using the ORBacus IDL to C++ translator, $IDL
        obv4=`echo $IDLVER | grep " 4"`
        obv4_1=`echo $IDLVER | grep " 4.1"`
        #echo obv4 = \"$obv4\"
        if test -n "$obv4"; then
          if test -n "$obv4_1"; then
            echo Using Orbacus 4.1
          else
            echo Using Orbacus 4.0
          fi
        else
          echo Unknown ORBabus version: quitting!
          exit
        fi
        CXXORB=ORBacus
        TESTINC=OB/CORBA.h
      else
        echo Unknown ORB: quitting!
        exit
      fi
    fi
  fi
fi

echo IDL = $IDL
if test -n "$IDLVER"; then
  if test -n "$isomni" -o -n "$istao"; then
    # This command works for tao and omniORB
    IDLVER=`$IDL -V 2>&1`
  else
    # This command is for orbacus and mico
    IDLVER=`$IDL --version 2>&1`
  fi
fi
echo IDLVER = $IDLVER



dnl ######################################################################
dnl
dnl Get the CORBA include files
dnl
dnl ######################################################################

AC_ARG_WITH(corba-incdir,
[  --with-corba-incdir=<corba-incdir-dir>  set location of the CORBA include directory],
CORBA_INCDIR="$withval")

dnl If not known, check in typical directories

if test -n "$CORBA_INCDIR"; then
  CORBA_INCPATH=$CORBA_INCDIR
else
  case $CXXORB in
    tao)
      CORBA_INCPATH=$TAO_ROOT/tao:$TAO_ROOT/include/tao
      CORBA_INCDIR2=$ACE_ROOT
      ;;

    omniORB)
      CORBA_INCPATH=/usr/local/omniORB4/include/omniORB4
      ;;
    MICO)
      MICO_ROOT=`dirname $IDL_DIR`
      CORBA_INCPATH=$MICO_ROOT/include
      AC_SUBST(MICO_VERSION)
      IDL="MICO_ROOT=$MICO_ROOT LD_LIBRARY_PATH=$MICO_ROOT/lib:$LD_LIBRARY_PATH $IDL"
      ;;
    ORBabus)
      CORBA_INCPATH=$HOME/$host/OB4-${COMPDIR}/include:$HOME/$host/OB4/include
      CORBA_INCPATH=$CORBA_INCPATH:/usr/local/OB4-${COMPDIR}/include:
      CORBA_INCPATH=$CORBA_INCPATH:/usr/local/OB4/include:
      CORBA_INCPATH=$CORBA_INCPATH:/loc/OB4-${COMPDIR}/include:/loc/OB4/include
      CORBA_INCPATH=$CORBA_INCPATH:/mfe/local/OB4-${COMPDIR}/include:/mfe/local/OB4/include
      CORBA_INCPATH=$CORBA_INCPATH:/local/OB4-${COMPDIR}/include:/local/OB4/include
      # We no longer support OB-3
      ;;
    *)
      # Well, we don't really know what's going on.  This should never happen.
      CORBA_INCPATH=`dirname $IDL_DIR`/include:$CORBA_INCPATH
      CORBA_INCDIR2=.
      ;;
  esac
fi

# It is dubious to apply this to all ORBs.

AC_SUBST(TAO_ROOT)
AC_SUBST(ACE_ROOT)
AC_SUBST(ACE_INCLUDE)
AC_SUBST(TAO_INCLUDE)
AC_SUBST(ACE_LIB_PATH)

AC_SUBST(MICO_ROOT)
AC_SUBST(CORBA_CXXFLAGS)

echo CORBA_INCPATH = $CORBA_INCPATH
TX_PATH_FILE(ABS_CORBA_INCDIR, $TESTINC,"", $CORBA_INCPATH)

if test -z "$ABS_CORBA_INCDIR"; then
  AC_MSG_ERROR(CORBA includes not found in $CORBA_INCPATH.  Use --with-corba-incdir to set the location of $TESINC.)
fi

CORBA_INCSUBDIR=`dirname $ABS_CORBA_INCDIR`
CORBA_INCDIR=`dirname $CORBA_INCSUBDIR`
if test -n "$ismico"; then
    CORBA_INCDIR=$CORBA_INCSUBDIR
fi
CORBA_DIR=`dirname $CORBA_INCDIR`

echo Corba include directory is $CORBA_INCDIR
echo CORBA_DIR=$CORBA_DIR

AC_SUBST(CORBA_INCDIR)
AC_SUBST(CORBA_INCDIR2)

dnl ######################################################################
dnl
dnl Get the JThreads++ include files if using Orbacus 4
dnl
dnl ######################################################################

#echo "######################################################################"
#echo "obv4 = $obv4"
#echo "######################################################################"

if test -z "$obv4"; then

  JTC_INCDIR=.

else

  AC_ARG_WITH(jtc-incdir,
  [  --with-jtc-incdir=<jthreads include directory>      to set location of JThreads/C++ headers],
  JTC_INCDIR="$withval")

  if test -n "$JTC_INCDIR"; then
    JTC_INCPATH=$JTC_INCDIR
  else
    JTC_INCPATH=$HOME/$host/jtc-${COMPDIR}/include:$HOME/$host/jtc/include
    JTC_INCPATH=$JTC_INCPATH:/usr/local/jtc-${COMPDIR}/include:
    JTC_INCPATH=$JTC_INCPATH:/usr/local/jtc/include:
    JTC_INCPATH=$JTC_INCPATH:/loc/jtc-${COMPDIR}/include:/loc/jtc/include
    JTC_INCPATH=$JTC_INCPATH:/mfe/local/jtc-${COMPDIR}/include:/mfe/local/jtc/include
    JTC_INCPATH=$JTC_INCPATH:/local/jtc-${COMPDIR}/include:/local/jtc/include
  fi

#  AC_PATH_PROGS(ABS_JTC_INCDIR, JTC/JTC.h,"", $JTC_INCPATH)
  TX_PATH_FILE(ABS_JTC_INCDIR, JTC/JTC.h,"", $JTC_INCPATH)
  if test -z "$ABS_JTC_INCDIR"; then
    AC_MSG_ERROR(JTC includes not found in $JTC_INCPATH.  Use --with-jtc-incdir to set the location of JTC/JTC.h.)
  fi

  JTC_INCDIR1=`dirname $ABS_JTC_INCDIR`
  JTC_INCDIR=`dirname $JTC_INCDIR1`
  JTC_DIR=`dirname $JTC_INCDIR`

fi

AC_SUBST(JTC_INCDIR)
AC_SUBST(JTC_DIR)

dnl ######################################################################
dnl
dnl Determine the correct CORBA library
dnl
dnl ######################################################################

AC_ARG_WITH(corba-libdir,
[  --with-corba-libdir=<corba library directory>      to set location of corba libraries],
CORBA_LIBDIR="$withval")

# Set in case not found or used

if test -n "$CORBA_LIBDIR"; then
  CORBA_LIBPATH=$CORBA_LIBDIR
else
  # Try to guess the corba library directory
  case $CXXORB in

    ORBacus)
      CORBA_LIBPATH=$CORBA_DIR/lib/${COMPDIR}:$CORBA_DIR/lib

      TX_PATH_FILE(ABS_CORBA_LIB, libOB.a libOB.so, "", $CORBA_LIBPATH)
      if test -z "$ABS_CORBA_LIB"; then
        AC_MSG_ERROR(Unable to find libOB.a or libOB.so in $CORBA_LIBPATH.  Set the appropriate directory using --with-corba-libdir)
      fi
      CORBA_LIBDIR=`dirname $ABS_CORBA_LIB`
      if test -n "$obv4"; then
        if test -n "$obv4_1"; then
          CORBA_LIB="-lOB"
        else
          CORBA_LIB="-lOB -lOBBiDir"
        fi
      else
        CORBA_LIB="-lOB"
      fi
      CORBA_LIBS="-L$CORBA_LIBDIR $CORBA_LIB"

      dnl Find the JThreads/C++ libraries
      if test -n "$obv4"; then
        AC_ARG_WITH(jtc-libdir,
          [  --with-jtc-libdir=<jthreads libdir>      to set location of JThreads/C++ libraries],
          JTC_LIBDIR="$withval")
        if test -n "$JTC_LIBDIR"; then
          JTC_LIBPATH=$JTC_LIBDIR
        else
          JTC_LIBPATH=$JTC_DIR/lib/${COMPDIR}:$JTC_DIR/lib
	# =$HOME/$host/jtc-${COMPDIR}/lib:$HOME/$host/jtc/lib/$COMPDIR:/usr/local/jtc-${COMPDIR}/lib:/usr/local/jtc/lib/$COMPDIR:/loc/jtc-${COMPDIR}/lib:/loc/jtc/lib/$COMPDIR:/local/jtc-${COMPDIR}/lib:/mfe/local/jtc-${COMPDIR}/lib:/mfe/local/jtc/lib/$COMPDIR
        fi
        TX_PATH_FILE(ABS_JTC_LIB, libJTC.a libJTC.so, "", $JTC_LIBPATH)
        if test -z "$ABS_JTC_LIB"; then
          AC_MSG_WARN(Unable to find libJTC.a or libJTC.so in $JTC_LIBPATH.  Set the appropriate directory using --with-jtc-libdir if needed)
        else
          JTC_LIBDIR=`dirname $ABS_JTC_LIB`
          JTC_LIB="-lJTC"
          JTC_LIBS="-L$JTC_LIBDIR $JTC_LIB"
          CORBA_LIB="$CORBA_LIB $JTC_LIB"
          CORBA_LIBS="$CORBA_LIBS $JTC_LIBS"
        fi
      fi
#   Minimal libraries for server and client apps.
      CORBA_SVR_LIBS="$CORBA_LIBS"
      CORBA_CLT_LIBS="$CORBA_LIBS"

      ;;

    MICO)
      CORBA_LIBPATH=$CORBA_DIR/lib:$HOME/$host/mico-${COMPDIR}/lib:$HOME/$host/mico/lib/$COMPDIR:/usr/local/mico-${COMPDIR}/lib:/usr/local/mico/lib/$COMPDIR:/loc/mico-${COMPDIR}/lib:/loc/mico/lib/$COMPDIR:/mfe/local/mico-${COMPDIR}/lib:/mfe/local/mico/lib/$COMPDIR

      TX_PATH_FILE(ABS_CORBA_LIB, libmico2.3.11.a libmicocoss2.3.11.so, "", $CORBA_LIBPATH)
      if test -z "$ABS_CORBA_LIB"; then
        AC_MSG_ERROR(Unable to find libmico.a or libmico.so in $CORBA_LIBPATH.  Set the appropriate directory using --with-corba-libdir)
      fi
      CORBA_LIBDIR=`dirname $ABS_CORBA_LIB`
      CORBA_LIB="-lmico2.3.11 -lmicocoss2.3.11"
      CORBA_LIBS="-L$CORBA_LIBDIR $CORBA_LIB"
#   Minimal libraries for server and client apps.
      CORBA_SVR_LIBS="$CORBA_LIBS"
      CORBA_CLT_LIBS="$CORBA_LIBS"
      ;;

    omniORB)
      CORBA_LIBPATH=/usr/local/omniORB4/lib/i586_linux_2.0_glibc2.1

      TX_PATH_FILE(ABS_CORBA_LIB, libomniORB4.a libomniORB4.so, "", $CORBA_LIBPATH)
      if test -z "$ABS_CORBA_LIB"; then
        AC_MSG_ERROR(Unable to find libomniORB4.a or libomniORB4.so in $CORBA_LIBPATH.  Set the appropriate directory using --with-corba-libdir)
      fi
      CORBA_LIBDIR=`dirname $ABS_CORBA_LIB`
      CORBA_LIB="-lomniORB4 -lbsd -lpthread -lomnithread -lomniDynamic"
      CORBA_LIBS="-L$CORBA_LIBDIR $CORBA_LIB"
#   Minimal libraries for server and client apps.
      CORBA_SVR_LIBS="$CORBA_LIBS"
      CORBA_CLT_LIBS="$CORBA_LIBS"
      ;;

    tao)
      CORBA_LIBPATH=$ACE_LIB_PATH
      TX_PATH_FILES(ABS_CORBA_LIB, libTAO.a libTAO.so libTAO.dylib, "", 
	$CORBA_LIBPATH)
      if test -z "$ABS_CORBA_LIB"; then
        AC_MSG_ERROR(Unable to find libTAO.a or libTAO.so or libTAO.dylib in $CORBA_LIBPATH.  Set the appropriate directory using --with-corba-libdir)
      fi
#    CORBA_LIBDIR=`dirname $ABS_CORBA_LIB`
      CORBA_LIBDIR=$ACE_LIB_PATH
      CORBA_LIB="-lTAO_CosNaming -lTAO_IFR_Client -lTAO_Svc_Utils -lTAO_IORTable -lTAO_BiDirGIOP -lTAO_DynamicInterface -lTAO_DynamicAny -lTAO_TypeCodeFactory -lTAO_DynamicAny -lTAO_Strategies -lTAO_PortableServer -lTAO_AnyTypeCode -lTAO -lACE"
#   Minimal libraries for server and client applications.
      CORBA_SVR_LIBS="-L$CORBA_LIBDIR -lTAO_PortableServer -lTAO_AnyTypeCode -lTAO -lACE -lTAO_Strategies"
      CORBA_CLT_LIBS="-L$CORBA_LIBDIR -lTAO -lACE"
      CORBA_LIBS="-L$CORBA_LIBDIR $CORBA_LIB"
      ;;

  esac
fi

dnl make substitutions
AC_SUBST(CORBA_LIBDIR)
AC_SUBST(CORBA_LIB)
AC_SUBST(CORBA_SVR_LIBS)
AC_SUBST(CORBA_CLT_LIBS)
AC_SUBST(CORBA_LIBS)
if test -z "$JTC_LIBDIR"; then JTC_LIBDIR=.; fi
AC_SUBST(JTC_LIBDIR)
AC_SUBST(JTC_LIB)

dnl ######################################################################
dnl
dnl    Write to summary file if defined
dnl
dnl ######################################################################
if test -n "$config_summary_file"; then
    echo                                         >> $config_summary_file
    echo "Using CORBA with"                      >> $config_summary_file
    echo "  CORBA_LIBDIR: $CORBA_LIBDIR"         >> $config_summary_file
    echo "  CORBA_LIB:   $CORBA_LIB"             >> $config_summary_file
    echo "  CORBA_INCDIR: $CORBA_INCDIR"         >> $config_summary_file
    echo "  CORBA_INCDIR2: $CORBA_INCDIR2"       >> $config_summary_file
    echo "  CORBA_SVR_LIBS:  $CORBA_SVR_LIBS"    >> $config_summary_file
    echo "  CORBA_CLT_LIBS:  $CORBA_CLT_LIBS"    >> $config_summary_file
    echo "  CORBA_LIBS:      $CORBA_LIBS"        >> $config_summary_file
    echo "  JTC_LIBDIR:      $JTC_LIBDIR"        >> $config_summary_file
    echo "  JTC_LIB:         $JTC_LIB"           >> $config_summary_file
fi
