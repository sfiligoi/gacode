dnl ######################################################################
dnl
dnl File:   fcflags.m4
dnl
dnl Purpose:
dnl         Given FC (f90) compiler, determine flags
dnl         See also fc.m4 and tx_fc_aux.m4
dnl         tx_fc_aux.m4 in particular defines the LDFLAGS for fortran
dnl         and handles determining the flag for including modules.
dnl
dnl Variables set here include
dnl      FC_LD
dnl      FC_DEFINE_FLAG to define for the preprocessor
dnl      FC_DBL_FLAG
dnl        Promote r4 to r8
dnl      FC_INT_FLAG
dnl        Promote i4 to i8
dnl      FC_WARN_FLAG
dnl        Turn on additional warnings. NEED WORK: defined only for lf95!
dnl      FC_OPT_FLAG
dnl        Normal optimization
dnl      FC_MACH_FLAG
dnl        Machine dependent flags for "ultra" optimization.  NEEDS WORK!!
dnl      FC_NOOPT_FLAG
dnl        o Sometimes better to explicitly turn off optimization
dnl      FC_DEBUG_FLAG
dnl        Include check array conformance if available
dnl      FC_OTHER_FLAG
dnl        Things should probably be used, but are not strictly necessary
dnl      FC_INCLUDE_FLAG
dnl        Because Absoft does not use the -I flag for include directories
dnl      FC_PIC_FLAG
dnl        Position-Independent Code (for shared libraries)
dnl      FCFLAGS
dnl      FC_VERSION
dnl      FC_LIBSUBDIR
dnl
dnl Description of differences between levels of optimization:
dnl   Standard: basic optimizations that are low-risk and do not lead to
dnl   very long compilation times or code bloat. Full compliance with the
dnl   IEEE-754 standard for floating point accuracy. No fancy tricks like
dnl   loop unrolling, inlining of non-inline functions, instruction reordering,
dnl   rewriting floating point expressions, or inter-procedural analysis. No
dnl   cpu-specific code. This generally corresponds to '-O2' or so.
dnl
dnl   Full: aggressive optimizations that might lead to longer compile
dnl   times and/or code size increase. Compliance with the IEEE standard
dnl   is relaxed. Function inlining permitted. Other fancy tricks also allowed
dnl   to some extend. No cpu-specific code. This generally corresponds to '-O3'.
dnl
dnl   Ultra: all caution to the wind. Compliance with IEEE goes out of the
dnl   window. Any trick the compiler provides may be used. Code can be
dnl   tailored to a specific cpu. Compiler options are highly platform-dependent
dnl
dnl NOTES:
dnl        o FC is the autotools awful variable name for any fortran
dnl          version beyond F77.  We use it interchangably with f90
dnl
dnl Version:      $Id: fcflags.m4 3784 2011-01-28 17:53:16Z pankin $
dnl
dnl Copyright 2001-2010, Tech-X Corporation.  Redistribution allowed provided
dnl this copyright statement remains intact.
dnl
dnl ######################################################################
dnl  consistency with flags.m4
dnl  some more work probably needed to ensure all of the flags are
dnl  accurate

# Change optimization specification.  Use with-optimization.
#
# OLD:
# Check for optimization. fullopt (-O2) is the default, so --enable-optimize
# is ignored, but must be listed here to be treated as a valid option by
# tx_check_args.m4.
#
# Definitions:
#
# --enable-minimalopt: -O
#
# --enable-optimize: The highest level of optimization such that the code
# remains valid for the processor family.
#
# --enable-fulloptimize: The highest level of optimization for the specific
# processor, rendering the code unusable on other processors in its family.
#
# --enable-ultraoptimize: Perform ALL optimizations, even ones which are not
# recommended and may possibly break the code.

dnl set the defaults
if test -z "$DOUBLE"; then DOUBLE=yes; fi
if test -z "$INT_DOUBLE"; then INT_DOUBLE=no; fi

if test -n "$FCFLAGS"; then SAVE_FCFLAGS="$FCFLAGS"; fi
if test -n "$FC_INCLUDE_FLAG"; then SAVE_FC_INCLUDE_FLAG="$FC_INCLUDE_FLAG"; fi
if test -n "$FC_DBL_FLAG"; then SAVE_FC_DBL_FLAG="$FC_DBL_FLAG"; fi
if test -n "$FC_OPT_FLAG"; then SAVE_FC_OPT_FLAG="$FC_OPT_FLAG"; fi
if test -n "$FC_MACH_FLAG"; then SAVE_FC_MACH_FLAG="$FC_MACH_FLAG"; fi
if test -n "$FC_FREE_FLAG"; then SAVE_FC_FREE_FLAG="$FC_FREE_FLAG"; fi
if test -n "$FC_FIXED_FLAG"; then SAVE_FC_FIXED_FLAG="$FC_FIXED_FLAG"; fi
if test -n "$FC_NOOPT_FLAG"; then SAVE_FC_NOOPT_FLAG="$FC_NOOPT_FLAG"; fi
if test -n "$FC_DEBUG_FLAG"; then SAVE_FC_DEBUG_FLAG="$FC_DEBUG_FLAG"; fi
if test -n "$FC_ALWAY_FLAG"; then SAVE_FC_ALWAY_FLAG="$FC_ALWAY_FLAG"; fi
if test -n "$FC_PIC_FLAG"; then SAVE_FC_PIC_FLAG="$FC_PIC_FLAG"; fi
if test -n "$FC_LDFLAGS"; then SAVE_FC_LDFLAGS="$FC_LDFLAGS"; fi

dnl ----------------------------------------------------------------------
dnl  allow flags or environment to overwrite variables
dnl ----------------------------------------------------------------------

AC_ARG_WITH(FCFLAGS,
  AC_HELP_STRING([--with-FCFLAGS=<desired fortran flags>],
  [Fortran flags]), FCFLAGS="$withval")

# New methodology for optimization, set string
AC_ARG_WITH(optimization,
  AC_HELP_STRING([--with-optimization=<opt level>],
    [Optimization level: one of debug, noopt, minimal, default, full, ultra]),
  OPTIMIZATION=${withval}, OPTIMIZATION=default)

AC_ARG_ENABLE(optimize,
  AC_HELP_STRING([--disable-optimize],[no longer valid]),
  AC_MSG_ERROR(--disable-optimize is no longer valid. Use --with-optimization=))

AC_ARG_ENABLE(debug,
  AC_HELP_STRING([--enable-debug],[no longer valid]),
  AC_MSG_ERROR(--enable-debug is no longer valid.  Use --with-optimization=))
if test "$DEBUG" = yes; then
  OPTIMIZATION="no"
fi

AC_ARG_ENABLE(fulloptimize,
  AC_HELP_STRING([--enable-fulloptimize], [no longer valid]),
  AC_MSG_ERROR(--enable-fulloptimize is no longer valid.  Use --with-optimization=))

AC_ARG_ENABLE(ultraoptimize,
  AC_HELP_STRING([--enable-ultraoptimize], [no longer valid]),
  AC_MSG_ERROR(--enable-ultraoptimize is no longer valid.  Use --with-optimization=))

# For backward compatibility
if test -z "$SERIALFC"; then
  SERIALFC=$SERIAL_F90
fi

echo "Setting the flags per system and FC compiler:" $SERIALFC

# Find real name of SERIALFC if a link.
case $SERIALFC in
  f95)
    TMP_SERIALFC=`which $SERIALFC`
    TMP_SERIALFC=`readlink -f $TMP_SERIALFC`
    TMP_SERIALFC=`basename $TMP_SERIALFC`
    if test -n "$TMP_SERIALFC"; then
      SERIALFC=$TMP_SERIALFC
    else
      AC_MSG_ERROR(fcflags.m4: unable to follow SERIALFC to its real name. Pls fix.)
    fi
    ;;
esac

dnl ----------------------------------------------------------------------
dnl Check on host and f90 since it is possible to have same compiler on
dnl different platforms with subtle differences
dnl ----------------------------------------------------------------------

processor=`uname -p`
echo SERIALFC = $SERIALFC
case $SERIALFC in

dnl Put xlf before lf as more restrictive
  xlf* | */xlf* | bgxlf* | */bgxlf*)
    dnl FC_DBL_FLAG="-qautodbl=dbl4"
    dnl FC_DBL_FLAG="-autodouble "
    dnl FC_DBL_FLAG="-qrealsize=8"
    # echo xlf found.
    FC_DBL_FLAG="-qrealsize=8 -qautodbl=dbl4"
    case "$sys" in
      AIX)
        FC_DEFINE_FLAG=-WF,-D
        ;;
      *)
dnl JRC: Not understood.  xlf_r on Intrepid does not this form in uestat?
        FC_DEFINE_FLAG=-WF,-D
        # FC_DEFINE_FLAG=-D
        ;;
    esac
    case $OPTIMIZATION in
      minimal) FC_OPT_FLAG="-O" ;;
      full)  FC_OPT_FLAG="-O3 -qstrict" ;;
      ultra) FC_OPT_FLAG="-O3 -qstrict" ;;
      *)   FC_OPT_FLAG="-O2" ;;
    esac
    FC_NOOPT_FLAG="-O0"
    FC_DEBUG_FLAG="-g"
# To find flush correctly
    FC_ALWAY_FLAG="-w -qextname=flush"
    case `uname` in
      AIX) # These don't work on BGP
        FC_ALWAY_FLAG="$FC_ALWAY_FLAG -qalign=natural"
        ;;
    esac
    FC_FIXED_FLAG="-qfixed=132"
    FC_FREE_FLAG="-qfree -qsuffix=f=f90"
    FC_PIC_FLAG="-qpic"
    FC_MACH_FLAG="-qarch=auto -qtune=auto -qcache=auto"
    FC_LIBSUBDIR=xlf
    ;;

  lf95* | */lf95*)
    # echo lf95 found.
    FC_DBL_FLAG="--dbl "
    FC_DEFINE_FLAG=-D
    case $OPTIMIZATION in
      full) FC_OPT_FLAG="-O3";;
      ultra) FC_OPT_FLAG="-O3";;
      *)  FC_OPT_FLAG="-O --prefetch 2" ;;
    esac
    FC_NOOPT_FLAG="-O0"
    FC_DEBUG_FLAG="-g --chk --trap --trace -DDEBUG"
    FC_ALWAY_FLAG="--ap --pca"
    FC_PIC_FLAG="-shared"
    FC_WARN_FLAG="--f95"
    dnl turn on Fortran 95 confirmation warning
    case $processor in
         athlon)  FC_MACH_FLAG="" ;;
         i686)  FC_MACH_FLAG="" ;;
    esac
    FC_VERSION=`$SERIALFC --version | head -1 | sed -e 's/.*Release //' -e 's/ .*//' -e 's/^L//' -e 's/[[^0-9.]].*//'`
    FC_VERSION_MM=`echo $FC_VERSION | sed 's/\.\([[0-9]]\).*/.\1/'`
    FC_LIBSUBDIR=lf95${FC_VERSION_MM}
    ;;

  *pgf90 | *pgf95)
    # echo pgf found.
    FC_DBL_FLAG="-r8 "
    FC_DEFINE_FLAG=-D
    dnl SEK: I think fast has all of these other optimizations but not sure
    dnl  FC_OPT_FLAG="-fast -Mcache_align -Munroll -Mdalign -Minline -Mvect=prefetch"
    case $OPTIMIZATION in
      full) FC_OPT_FLAG="-fast";;
      ultra) FC_OPT_FLAG="-fast";;
      *)  FC_OPT_FLAG="-O2";;
    esac
    FC_NOOPT_FLAG=""
    FC_DEBUG_FLAG="-g -C -Mbounds -Minfo -DDEBUG"
    dnl Use little endian for binary compatibility
    FC_ALWAY_FLAG="-byteswapio -Mextend"
    case $processor in
         athlon)  FC_MACH_FLAG="" ;;
         i686)  FC_MACH_FLAG="" ;;
    esac
    FC_LIBSUBDIR=pgf
    ;;

  *pathf90* | *pathf95*)
    AC_DEFINE(NoSystem,[],"Define whether the system call is supported")
    FC_DBL_FLAG="-r8 "
    FC_DEFINE_FLAG=-D
    dnl O3 seems to cause problems with NIMROD
    case $OPTIMIZATION in
      full) FC_OPT_FLAG="-O2";;
      ultra) FC_OPT_FLAG="-O2";;
      *)  FC_OPT_FLAG="-O2";;
    esac
    FC_NOOPT_FLAG=""
    FC_DEBUG_FLAG="-g -DDEBUG"
    dnl the -C is for bounds checking.  Latest version on pathf90
    dnl  doesn't work well
    dnl FC_DEBUG_FLAG="-g -C -DDEBUG"
    FC_ALWAY_FLAG="-extend-source -byteswapio -fno-second-underscore"
    FC_LDFLAGS="-Wl,--warn-unresolved-symbols"
    case $processor in
         athlon) FC_MACH_FLAG="" ;;
         x86_64) FC_MACH_FLAG="-march=opteron -mtune=opteron -msse3" ;;
    esac
    FC_LIBSUBDIR=pathf
    FC_PIC_FLAG=" "
    ;;
  *ifort)
    FC_DBL_FLAG="-autodouble "
    FC_DEFINE_FLAG=-D
    case $OPTIMIZATION in
      full) FC_OPT_FLAG="-O3";;
      ultra) FC_OPT_FLAG="-O3 -fno-alias";;
      *)  FC_OPT_FLAG="-O2";;
    esac
    FC_NOOPT_FLAG="-O0"
    FC_DEBUG_FLAG="-g -CB -debug inline_debug_info -DDEBUG"
    FC_ALWAY_FLAG="-fltconsistency -convert big_endian -132"
       dnl FC_ALWAY_FLAG="-cm -w -w95"
    FC_PIC_FLAG="-fpic"
    ##FC_MACH_FLAG="-O2 -tpp7 -xW"
    ## Assuming Core processors
    FC_MACH_FLAG="-axP"
    ;;
  *g95)
    FC_DBL_FLAG=""
    FC_DEFINE_FLAG=-D
    case $OPTIMIZATION in
      full) FC_OPT_FLAG="-O3" ;;
      ultra) FC_OPT_FLAG="-O5" ;;
      *)  FC_OPT_FLAG="-O2" ;;
    esac
    FC_NOOPT_FLAG="-O0"
    FC_DEBUG_FLAG="-g -DDEBUG"
    FC_ALWAY_FLAG="-fno-second-underscore"
    FC_MACH_FLAG=""
    ;;
  *gfortran*)
    FC_DBL_FLAG="-fdefault-real-8 -fdefault-double-8"
    FC_DEFINE_FLAG=-D
    FC_INT_FLAG=""
    dnl FC_INT_FLAG="-fdefault-integer-8"
    case $OPTIMIZATION in
      full)  FC_OPT_FLAG="-O3" ;;
      ultra) FC_OPT_FLAG="-O3" ;;
      *)     FC_OPT_FLAG="-O2 -g" ;;
    esac
    FC_NOOPT_FLAG="-O0"
    FC_DEBUG_FLAG="-g -fbounds-check -DDEBUG"
    dnl FC_ALWAY_FLAG="-fno-second-underscore"
    dnl FC_ALWAY_FLAG="-ffixed-line-length-132 -fimplicit-none -Wno-tabs"
dnl no-tabs breaks on Intrepid
dnl the fixed line length breaks nimrod
    if test "$USE_NIMFLAGS" = true; then
    FC_ALWAY_FLAG="-fno-second-underscore -fconvert=big-endian"
    else
    FC_ALWAY_FLAG="-ffixed-line-length-132"
    fi
    FC_MACH_FLAG=""
# Do not use SERIALFC, as it may not be in one's path
    FC_VERSION=`$FC --version | sed -n '1p' | sed -e 's/^.*GCC. //' -e 's/ .*$//'`
# Double brackets as shell swallows outer set
    FC_VERSION_MM=`echo $FC_VERSION | sed 's/\.[[0-9]]*$//'`
    FC_LIBSUBDIR=gfortran${FC_VERSION_MM}
    ;;
  *)
    FC_DBL_FLAG=""
    FC_DEFINE_FLAG=-D
    case $OPTIMIZATION in
      full) FC_OPT_FLAG="-O" ;;
      ultra) FC_OPT_FLAG="-O" ;;
      *)  FC_OPT_FLAG="-O" ;;
    esac
    FC_NOOPT_FLAG=""
    FC_DEBUG_FLAG="-g -DDEBUG"
    FC_ALWAY_FLAG=""
    FC_MACH_FLAG=""
    ;;
esac

dnl ----------------------------------------------------------------------
dnl Set up default flags by setting variables that haven't been defined
dnl ----------------------------------------------------------------------
if test -z "$FC_INCLUDE_FLAG"; then FC_INCLUDE_FLAG="-I"; fi

dnl Can't test on all platforms, so this could cause errors.
if test -z "$FC_PIC_FLAG"; then FC_PIC_FLAG="-fpic"; fi

dnl ----------------------------------------------------------------------
dnl  Allow overriding
dnl ----------------------------------------------------------------------
if test -n "$SAVE_FC_INCLUDE_FLAG"; then FC_INCLUDE_FLAG="$SAVE_FC_INCLUDE_FLAG"; fi
if test -n "$SAVE_FC_DBL_FLAG"; then FC_DBL_FLAG="$SAVE_FC_DBL_FLAG"; fi
if test -n "$SAVE_FC_OPT_FLAG"; then FC_OPT_FLAG="$SAVE_FC_OPT_FLAG"; fi
if test -n "$SAVE_FC_MACH_FLAG"; then FC_MACH_FLAG="$SAVE_FC_MACH_FLAG"; fi
if test -n "$SAVE_FC_FREE_FLAG"; then FC_FREE_FLAG="$SAVE_FC_FREE_FLAG"; fi
if test -n "$SAVE_FC_FIXED_FLAG"; then FC_FIXED_FLAG="$SAVE_FC_FIXED_FLAG"; fi
if test -n "$SAVE_FC_NOOPT_FLAG"; then FC_NOOPT_FLAG="$SAVE_FC_NOOPT_FLAG"; fi
if test -n "$SAVE_FC_DEBUG_FLAG"; then FC_DEBUG_FLAG="$SAVE_FC_DEBUG_FLAG"; fi
if test -n "$SAVE_FC_ALWAY_FLAG"; then FC_ALWAY_FLAG="$SAVE_FC_ALWAY_FLAG"; fi
if test -n "$SAVE_FC_PIC_FLAG"; then FC_PIC_FLAG="$SAVE_FC_PIC_FLAG"; fi
if test -n "$SAVE_FC_LDFLAGS"; then FC_LDFLAGS="$SAVE_FC_LDFLAGS"; fi

dnl ----------------------------------------------------------------------
dnl   Set up fc
dnl ----------------------------------------------------------------------

if test -n "$FC"; then
  if test -z "$FC_LD"; then
      FC_LD=$FC
  fi
fi
if test -z "$FC_LDFLAGS"; then
    FC_LDFLAGS=""
fi

dnl ----------------------------------------------------------------------
dnl   Set flags based on input.
dnl ----------------------------------------------------------------------

FCFLAGS="$FC_ALWAY_FLAG $FC_PIC_FLAG"

if test "$DOUBLE" = yes; then
   FCFLAGS="$FCFLAGS $FC_DBL_FLAG"
fi
if test "$INT_DOUBLE" = yes; then
   FCFLAGS="$FCFLAGS $FC_INT_FLAG"
fi

case $OPTIMIZATION in
  full | ultra)
    AC_MSG_WARN(May generate code for this processor only.)
    FCFLAGS="$FCFLAGS $FC_OPT_FLAG $FC_MACH_FLAG"
    ;;
  noopt) FCFLAGS="$FCFLAGS $FC_NOOPT_FLAG";;
  debug) FCFLAGS="$FCFLAGS $FC_DEBUG_FLAG";;
  *) FCFLAGS="$FCFLAGS $FC_OPT_FLAG";;
esac

dnl ----------------------------------------------------------------------
dnl FCFLAGS at the command-line is additive
dnl ----------------------------------------------------------------------
if test -n "$SAVE_FCFLAGS"; then FCFLAGS="$FCFLAGS $SAVE_FCFLAGS"; fi

dnl ----------------------------------------------------------------------
dnl  AC_SUBST everything to allow fine-grained control of compilation
dnl ----------------------------------------------------------------------

AC_SUBST(FCFLAGS)
AC_SUBST(FC_ALWAY_FLAG)
AC_SUBST(FC_DBL_FLAG)
AC_SUBST(FC_DEBUG_FLAG)
AC_SUBST(FC_DEFINE_FLAG)
AC_SUBST(FC_FIXED_FLAG)
AC_SUBST(FC_FREE_FLAG)
AC_SUBST(FC_INCLUDE_FLAG)
AC_SUBST(FC_LD)
AC_SUBST(FC_LDFLAGS)
AC_SUBST(FC_LIBSUBDIR)
AC_SUBST(FC_MACH_FLAG)
AC_SUBST(FC_NOOPT_FLAG)
AC_SUBST(FC_OPT_FLAG)
AC_SUBST(FC_PIC_FLAG)
AC_SUBST(FC_VERSION)
AC_SUBST(FC_WARN_FLAG)
AC_DEFINE_UNQUOTED([FC_VERSION], "$FC_VERSION", "Fortran compiler version")

dnl ----------------------------------------------------------------------
dnl  Print to config summary file if defined
dnl ----------------------------------------------------------------------
if test -n "$config_summary_file"; then
  echo                                        >> $config_summary_file
  echo "Fortran compiler and flags:"          >> $config_summary_file
  TX_PRINT_VAR(FC)
  TX_PRINT_VAR(FCFLAGS)
  TX_PRINT_VAR(FC_VERSION)
  TX_PRINT_VAR(FC_LIBSUBDIR)
fi

