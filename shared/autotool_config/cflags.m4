dnl ######################################################################
dnl
dnl File:	flags.m4
dnl
dnl Purpose:	To determine the flags for C
dnl
dnl SEK: This is done completely differently than Tech-X config
dnl
dnl Version:	$Id: cflags.m4 3593 2010-07-13 20:14:57Z kruger $
dnl
dnl ######################################################################

dnl ----------------------------------------------------------------------
dnl  allow flags or environment to overwrite variables
dnl ----------------------------------------------------------------------
AC_ARG_WITH(CFLAGS,
    [  --with-CFLAGS="MyCFLAGS"   to set C flags with CC],
    CFLAGS="$withval")

if test -n "$CFLAGS"; then
   dnl AC_PROG_CC can define CFLAGS
   if test "$CFLAGS" != "-g"; then
      SAVE_CFLAGS="$CFLAGS"
   fi
fi

dnl ----------------------------------------------------------------------
dnl Determine flags based on host and CC compiler
dnl ----------------------------------------------------------------------
case "$host" in
 *-linux-gnu)
  processor=`uname -p`
  case $CC in
     *gcc*) CFLAGS=""
            CFLAGS_DBG="-g"
            CFLAGS_OPT=""
	;;
     *pgcc) CFLAGS=""
            CFLAGS_DBG="-g"
            CFLAGS_OPT=""
	;;
     *pathcc*) CFLAGS="-DAdd__ -fno-unsafe-math-optimizations"
            CFLAGS_DBG="-g"
  	      case $processor in
	         athlon) CFLAGS_OPT="-lm -O2 -march=opteron -mcpu=opteron -msse3" ;;
	         i686)  CFLAGS_OPT="" ;;
	      esac
	;;
     *icc) CFLAGS=""
            CFLAGS_DBG="-g"
            CFLAGS_OPT=""
	;;
     *) CFLAGS=""
        CFLAGS_DBG="-g"
        CFLAGS_OPT=""
	;;
  esac
  ;;

 *apple-darwin*) CFLAGS=""; CFLAGS_DBG="-g"; CFLAGS_OPT="" ;;

 *-sgi-irix*)            CFLAGS=""; CFLAGS_DBG="-g"; CFLAGS_OPT="" ;;
 hppa*-hp-hpux*)      CFLAGS=""; CFLAGS_DBG="-g"; CFLAGS_OPT="" ;;
 alpha*-cray-unicos*) CFLAGS=""; CFLAGS_DBG="-g"; CFLAGS_OPT="" ;;

 *-*-aix*)
     case "$host" in
       powerpc64-*)
	   CFLAGS="-q64 -DALL_SOURCE -DNoChange -qlanglvl=ansi"
         CFLAGS_DBG="-g";
	   CFLAGS_OPT="-D__RS6000" ;;
       *)
	   CFLAGS="-DALL_SOURCE -DNoChange -qlanglvl=ansi"
         CFLAGS_DBG="-g"; CFLAGS_OPT="" ;;
	esac
	;;
esac

dnl ----------------------------------------------------------------------
dnl   Set flags based on input.  Debug flag overwrites optimization flag
dnl ----------------------------------------------------------------------

if test -n "$CC"; then
     case $OPTIMIZED in
       yes) CFLAGS="$CFLAGS $CFLAGS_OPT";;
       ultra) CFLAGS="$CFLAGS $CFLAGS_OPT";;
     esac

     if test "$DEBUG" = yes; then
        CFLAGS="$CFLAGS $CFLAGS_DBG"
     fi
fi

dnl Command line overwriate everything
if test -n "$SAVE_CFLAGS"; then CFLAGS="$SAVE_CFLAGS"; fi
dnl ----------------------------------------------------------------------
dnl  AC_SUBST everything to allow fine-grained control of compilation
dnl ----------------------------------------------------------------------
AC_SUBST(CFLAGS)
AC_SUBST(CFLAGS_DBG)
AC_SUBST(CFLAGS_OPT)

