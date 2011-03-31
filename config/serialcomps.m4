dnl ######################################################################
dnl
dnl File:	serialcomps.m4
dnl
dnl Purpose:	Determine the serial compilers for C, C++, Fortran.
dnl             Absolute path, e.g., ABS_SERIALCXX, and basename, SERIALCXX
dnl		to be found.
dnl		Suppose to find the unique name (e.g., sunCC instead of CC)
dnl		from any compiler wrapper or link.
dnl
dnl Version:	$Id: serialcomps.m4 3669 2010-09-09 22:11:27Z kruger $
dnl
dnl Copyright 2001-2010, Tech-X Corporation.  Redistribution allowed provided
dnl this copyright statement remains intact.
dnl
dnl ######################################################################

echo "Determining the serial compilers."

dnl ######################################################################
dnl
dnl Determine the serial C++ compiler
dnl
dnl ######################################################################

if test -n "$CXX"; then

  CXXBASE=`echo $CXX | sed 's/ .*$//'`
  CXXBASE=`basename $CXXBASE`
  case $CXXBASE in
    CC)

# Assume existence of PE_ENV means CRAY
      case "$PE_ENV" in
        PATHSCALE)
          if test -z "$SERIALCXX"; then
            SERIALCXX=pathCC
          fi
          ;;
        PGI)
          if test -z "$SERIALCXX"; then
            SERIALCXX=pgCC
          fi
          ;;
        GNU)
          if test -z "$SERIALCXX"; then
            SERIALCXX=g++
          fi
          ;;
      esac
      if test -n "$SERIALCXX"; then
        ABS_SERIALCXX=`which $CXX`
      else
# Otherwise check for sun
        isSun=`$CXX -V 2>&1 | sed -n '1p' | grep " Sun "`
        if test -n "$isSun"; then
          SERIALCXX=sunCC
          ABS_SERIALCXX=`which $CXX`
          ABS_SERIALCXX=`dirname $ABS_SERIALCXX`/sunCC
        else
          AC_MSG_ERROR(Unable to determine the real compiler inside of CC)
        fi
      fi

      if false; then
      knowCC=false
      CC -show 1>/dev/null 2>&1
# Really pgCC compiler?
      if test $? = 0; then
        ispgi=`CC -show 2>/dev/null | grep PGILD`
        if test -n "$ispgi"; then
          SERIALCXX=pgCC
          ABS_SERIALCXX=`CC -show 2>/dev/null | grep DRIVERDIR | sed 's/^.*=//'`
          ABS_SERIALCXX=`echo $ABS_SERIALCXX | sed -e 's/^ *//' -e 's/  *$//'`
          ABS_SERIALCXX=${ABS_SERIALCXX}/pgCC
          knowCC=true
        fi
# Really gcc?
        if ! $knowCC; then
          isGcc=`$CXX --version 2>&1 | grep "g++"`
          if test -n "$isGcc"; then
            SERIALCXX=g++
            ABS_SERIALCXX=`which $CXX`
            ABS_SERIALCXX=`dirname $ABS_SERIALCXX`/g++
            knowCC=true
          fi
        fi
# Really path?
        if ! $knowCC; then
          isPath=`$CXX --version 2>&1 | grep "PathScale"`
          if test -n "$isPath"; then
            SERIALCXX=pathCC
            ABS_SERIALCXX=`which $CXX`
            ABS_SERIALCXX=`dirname $ABS_SERIALCXX`/CC
            knowCC=true
          fi
        fi
      fi
# Really sunCC?
      if ! $knowCC; then
        isSun=`$CXX -V 2>&1 | sed -n '1p' | grep " Sun "`
        if test -n "$isSun"; then
          SERIALCXX=sunCC
          ABS_SERIALCXX=`which $CXX`
          ABS_SERIALCXX=`dirname $ABS_SERIALCXX`/sunCC
          knowCC=true
        fi
      fi
      fi

      if ! $knowCC; then
        AC_MSG_ERROR(Unable to determine the real compiler inside of CC)
      fi
      ;;

    mpiCC | mpicxx | mpic++)
      MPI_OPTS=`$CXX -show`
      ABS_SERIALCXX=`echo $MPI_OPTS | sed 's/ .*$//'`
      dnl echo SERIALCXX = $SERIALCXX
      if test -z "$ABS_SERIALCXX"; then
        AC_MSG_ERROR(Unable to determine the real compiler inside of $CXX)
      fi
      SERIALCXX=`basename $ABS_SERIALCXX`
      ;;
    mpCC)
      SERIALCXX=xlC
      ABS_SERIALCXX=`which $SERIALCXX`
      ;;
    mpCC_r | mpixlcxx_r)
      SERIALCXX=xlC_r
      ABS_SERIALCXX=`which $SERIALCXX`
      ;;
    mpiicpc)
      SERIALCXX=icpc
      ABS_SERIALCXX=`which $SERIALCXX`
      ;;
    mpipathCC)
      SERIALCXX=pathCC
      ABS_SERIALCXX=`which $SERIALCXX`
      ;;
    mpipgCC)
      SERIALCXX=pgCC
      ABS_SERIALCXX=`which $SERIALCXX`
      ;;
    *)
      SERIALCXX=$CXXBASE
      ABS_SERIALCXX=`echo $CXX | sed 's/ .*$//'`
      ABS_SERIALCXX=`which $ABS_SERIALCXX`
      ;;
  esac
fi

AC_SUBST(CXXBASE)
AC_SUBST(SERIALCXX)
AC_SUBST(ABS_SERIALCXX)
echo Serial C++ compiler is \`$SERIALCXX\'
echo Absolute path to C++ compiler is \`$ABS_SERIALCXX\'

dnl ######################################################################
dnl
dnl Determine the serial C compiler
dnl
dnl ######################################################################

dnl No conditional, as we always determine the C compiler

# In case flags included
CCBASE=`echo $CC | sed 's/ .*$//'`
CCBASE=`basename $CCBASE`
case $CCBASE in

  cc)

# Assume existence of PE_ENV means CRAY
    case "$PE_ENV" in
      PATHSCALE)
        if test -z "$SERIALCC"; then
          SERIALCC=pathcc
        fi
        ;;
      PGI)
        if test -z "$SERIALCC"; then
          SERIALCC=pgcc
        fi
        ;;
      GNU)
        if test -z "$SERIALCC"; then
          SERIALCC=gcc
        fi
        ;;
    esac
    if test -n "$SERIALCC"; then
      ABS_SERIALCC=`which $CC`
    else
# Otherwise check for sun
      isSun=`$CC -V 2>&1 | sed -n '1p' | grep " Sun "`
      if test -n "$isSun"; then
        SERIALCC=suncc
        ABS_SERIALCC=`which $CC`
        ABS_SERIALCC=`dirname $ABS_SERIALCC`/suncc
      else
        AC_MSG_ERROR(Unable to determine the real compiler inside of cc)
      fi
    fi

    if false; then
# sh does not have pipestatus
    $CC -show 1>/dev/null 2>&1
# Really pgcc?
    knowcc=false
    if test $? = 0; then
      ispgi=`cc -show 2>/dev/null | grep PGILD`
      if test -n "$ispgi"; then
        SERIALCC=pgcc
        ABS_SERIALCC=`CC -show 2>/dev/null | grep DRIVERDIR | sed 's/^.*=//'`
        ABS_SERIALCC=`echo $ABS_SERIALCC | sed -e 's/^ *//' -e 's/  *$//'`
        ABS_SERIALCC=${ABS_SERIALCC}/pgcc
        knowcc=true
      fi
# Really gcc?
      if ! $knowcc; then
        isGcc=`$CC -V 2>&1 | grep "gcc"`
        if test -n "$isGcc"; then
          SERIALCC=gcc
          ABS_SERIALCC=`which $CC`
          ABS_SERIALCC=`dirname $ABS_SERIALCC`/gcc
          knowcc=true
        fi
      fi
      if ! $knowcc; then
        isPath=`$CC --version 2>&1 | grep "PathScale"`
        if test -n "$isPath"; then
          SERIALCC=pathcc
          ABS_SERIALCC=`which $CC`
          ABS_SERIALCC=`dirname $ABS_SERIALCC`/cc
          knowcc=true
        fi
      fi
    fi
# Really suncc?
    if ! $knowcc; then
      isSun=`$CC -V 2>&1 | sed -n '1p' | grep " Sun "`
      if test -n "$isSun"; then
        SERIALCC=suncc
        ABS_SERIALCC=`which $CC`
        ABS_SERIALCC=`dirname $ABS_SERIALCC`/suncc
        knowcc=true
      fi
    fi
# Give up
    if ! $knowcc; then
      AC_MSG_ERROR(Unable to determine the real compiler inside of cc)
    fi
    fi

    ;;
  mpicc)
    MPI_OPTS=`$CC -show`
    ABS_SERIALCC=`echo $MPI_OPTS | sed 's/ .*$//'`
    dnl echo SERIALCC = $SERIALCC
    if test -z "$ABS_SERIALCC"; then
      AC_MSG_ERROR(Unable to determine the real compiler inside of $CC)
    fi
    SERIALCC=`basename $ABS_SERIALCC`
    ;;
  mpcc)
    SERIALCC=xlc
    ABS_SERIALCC=`which $SERIALCC`
    ;;
  mpcc_r | mpixlc_r)
    SERIALCC=xlc_r
    ABS_SERIALCC=`which $SERIALCC`
    ;;
  mpiicc)
    SERIALCC=icc
    ABS_SERIALCC=`which $SERIALCC`
    ;;
  mpipathcc)
    SERIALCC=pathcc
    ABS_SERIALCC=`which $SERIALCC`
    ;;
  mpipgcc)
    SERIALCC=pgcc
    ABS_SERIALCC=`which $SERIALCC`
    ;;
  *)
    SERIALCC=$CCBASE
    ABS_SERIALCC=`echo $CC | sed 's/ .*$//'`
    ABS_SERIALCC=`which $ABS_SERIALCC`
    ;;
esac

AC_SUBST(CCBASE)
AC_SUBST(SERIALCC)
AC_SUBST(ABS_SERIALCC)
echo Serial C compiler is \`$SERIALCC\'
echo Absolute path to C compiler is \`$ABS_SERIALCC\'

dnl ######################################################################
dnl
dnl JRC 4Aug09: Fortran compilers need to be made compliant with above!
dnl NEED ABS_SERIALFC, ABS_SERIALF77
dnl
dnl ######################################################################

dnl ######################################################################
dnl
dnl Determine the serial Fortran compiler
dnl
dnl ######################################################################

if test -z "$FC"; then
  NO_FORTRAN=true
fi
dnl SEK -- serialcomps.m4 is always going to be called before
dnl tx_fc_aux.m4 so this will be overwitten if needed; however,
dnl by removing the if statement we remove the 
dnl "AM_CONDITIONAL in if-branch" problem.
dnl if test "$NO_FORTRAN" = true; then
  AM_CONDITIONAL(TX_FORTRAN_MODNAMEISCAP, false)
dnl fi

if test "$NO_FORTRAN" != true; then
  echo Free-form Fortran compiler is \`$FC\'
  if test -n "$FC"; then
    FCBASE=`echo $FC | sed 's/ .*$//'`
    FCBASE=`basename $FCBASE`
    case $FCBASE in
      ftn)

# Assume existence of PE_ENV means CRAY
        case "$PE_ENV" in
          PATHSCALE)
            if test -z "$SERIALFC"; then
              SERIALFC=pathf90
            fi
            ;;
          PGI)
            if test -z "$SERIALFC"; then
              SERIALFC=pgf90
            fi
            ;;
          GNU)
            if test -z "$SERIALFC"; then
              SERIALFC=gfortran
            fi
            ;;
        esac
        if test -n "$SERIALFC"; then
          ABS_SERIALFC=`which $FC`
        else
# Otherwise check for sun
          isSun=`$FC -V 2>&1 | sed -n '1p' | grep " Sun "`
          if test -n "$isSun"; then
# This may be wrong.  Will fix when have a machine.
            SERIALFC=sunftn
            ABS_SERIALFC=`which $FC`
            ABS_SERIALFC=`dirname $ABS_SERIALFC`/sunftn
          else
            AC_MSG_ERROR(Unable to determine the real compiler inside of ftn)
          fi
        fi

        if false; then
# sh does not have pipestatus
        knowftn=false
        ftn -show 1>/dev/null 2>&1
        if test $? = 0; then
          ispgi=`ftn -show 2>/dev/null | grep PGILD`
          if test -n "$ispgi"; then
            SERIALFC=pgf90
            knowftn=true
          fi
          # echo After pgi, knowftn = $knowftn
          if ! $knowftn; then
            isGcc=`$FC -V 2>&1 | grep "gfortran"`
            # echo isGcc = $isGcc
            if test -n "$isGcc"; then
              SERIALFC=gfortran
              knowftn=true
            fi
          fi
          # echo After gcc, knowftn = $knowftn
          if ! $knowftn; then
            isPath=`$FC -V 2>&1 | grep "pathf90"`
            # echo isPath = $isPath
            if test -n "$isPath"; then
              SERIALFC=pathf90
              knowftn=true
            fi
          fi
          if ! $knowftn; then
            AC_MSG_ERROR(Unable to determine the real compiler inside of ftn)
          fi
        else
          SERIALFC=ftn
        fi
        fi

        ;;
      mpif9? | mpif03)
        MPI_OPTS=`$FC -show`
        SERIALFC=`echo $MPI_OPTS | sed 's/ .*$//'`
        if test -z "$SERIALFC"; then
          AC_MSG_ERROR(Unable to determine the real compiler inside of $F90)
        fi
        SERIALFC=`basename $SERIALFC`
        ;;
      mpipgf90)
        SERIALFC=pgf90
        ;;
      mpxlf | mpxlf9? | mpxlf03)
        SERIALFC=xlf
        ;;
      mpxlf_r | mpxlf9?_r | mpxlf03_r | mpixlf_r | mpixlf9?_r | mpixlf03_r)
        SERIALFC=xlf_r
        ;;
      *)
        SERIALFC=$FCBASE
        ;;
    esac
  fi
  echo Serial free-form Fortran compiler \(FC\) is \`$SERIALFC\'
fi
AC_SUBST(FCBASE)
AC_SUBST(SERIALFC)

dnl ######################################################################
dnl
dnl Determine the serial Fortran compiler
dnl
dnl ######################################################################

if test "$NO_FORTRAN" != true; then
  echo Fixed-form Fortran compiler is \`$F77\'
  if test -n "$F77"; then
    F77BASE=`echo $F77 | sed 's/ .*$//'`
    F77BASE=`basename $F77BASE`
    case $F77BASE in
      ftn)

# Assume existence of PE_ENV means CRAY
        case "$PE_ENV" in
          PATHSCALE)
            if test -z "$SERIALF77"; then
              SERIALF77=pathf90
            fi
            ;;
          PGI)
            if test -z "$SERIALF77"; then
              SERIALF77=pgf77
            fi
            ;;
          GNU)
            if test -z "$SERIALF77"; then
              SERIALF77=gfortran
            fi
            ;;
        esac
        if test -n "$SERIALF77"; then
          ABS_SERIALF77=`which $F77`
        else
# Otherwise check for sun
          isSun=`$F77 -V 2>&1 | sed -n '1p' | grep " Sun "`
          if test -n "$isSun"; then
# This may be wrong.  Will fix when have a machine.
            SERIALF77=sunftn
            ABS_SERIALF77=`which $F77`
            ABS_SERIALF77=`dirname $ABS_SERIALF77`/sunftn
          else
            AC_MSG_ERROR(Unable to determine the real compiler inside of ftn)
          fi
        fi

        if false; then
# sh does not have pipestatus
        ftn -show 1>/dev/null 2>&1
        if test $? = 0; then
          ispgi=`$F77 -show 2>/dev/null | grep PGILD`
          if test -n "$ispgi"; then
            SERIALF77=pgf90
            knowftn=true
          fi
          # echo After pgi, knowftn = $knowftn
          if ! $knowftn; then
            isGcc=`$F77 -V 2>&1 | grep "gfortran"`
            # echo isGcc = $isGcc
            if test -n "$isGcc"; then
              SERIALF77=gfortran
              knowftn=true
            fi
          fi
          # echo After gcc, knowftn = $knowftn
          if ! $knowftn; then
            isPath=`$F77 -V 2>&1 | grep "pathf90"`
            # echo isPath = $isPath
            if test -n "$isPath"; then
              SERIALF77=pathf90
              knowftn=true
            fi
          fi
          if ! $knowftn; then
            AC_MSG_ERROR(Unable to determine the real compiler inside of ftn)
          fi
        else
          SERIALF77=ftn
        fi
        fi

        ;;
      mpif*)
        MPI_OPTS=`$F77 -show`
        SERIALF77=`echo $MPI_OPTS | sed 's/ .*$//'`
        if test -z "$SERIALF77"; then
          AC_MSG_ERROR(Unable to determine the real compiler inside of $F77)
        fi
        SERIALF77=`basename $SERIALF77`
        ;;
      mpipgf*?)
        SERIALFC=pgf77
        ;;
dnl JRC: For xl, F77 should never be anything but the xlf or xlf_r compiler
      mpxlf77 | mpixlf77)
dnl SriV        SERIALF77=xlf77
        SERIALF77=xlf
        ;;
      mpxlf77_r | mpixlf77_r)
dnl Sriv        SERIALF77=xlf77_r
        SERIALF77=xlf_r
        ;;
      *)
        SERIALF77=$F77BASE
        ;;
    esac
  fi
  echo Serial free-form Fortran compiler \(F77\) is \`$SERIALF77\'
fi
dnl echo F77 = $F77, SERIALF77 = $SERIALF77
AC_SUBST(F77BASE)
AC_SUBST(SERIALF77)

dnl ######################################################################
dnl
dnl Write to summary file
dnl
dnl ######################################################################

echo >>$config_summary_file
echo Serial compilers >>$config_summary_file
TX_PRINT_VAR(SERIALCC)
TX_PRINT_VAR(SERIALCXX)
# echo "  SERIALCC:       $SERIALCC"        >> $config_summary_file
# echo "  SERIALCXX:      $SERIALCXX"       >> $config_summary_file
if test "$NO_FORTRAN" != true; then
  TX_PRINT_VAR(SERIALFC)
  TX_PRINT_VAR(SERIALF77)
  # echo "  SERIALFC:       $SERIALFC"        >> $config_summary_file
  # echo "  SERIALF77:      $SERIALF77"       >> $config_summary_file
fi

