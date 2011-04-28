dnl ######################################################################
dnl
dnl File:     tx_petsc3.m4
dnl
dnl Purpose:  Determine if PETSc is installed.  This is a remake
dnl           Tried Satish method, but it just gave too much crap -
dnl           unneeded directories and such.  So will now just do
dnl           the traditional..
dnl
dnl Variables defined in AC_SUBST
dnl   Basics:			Include Directories:
dnl     PETSC_DIR		  PETSC_INCDIR
dnl     PETSC_ARCH		  PETSC_INCDIR_ARCH
dnl   Library:			Library Info:
dnl     PETSC_LIB		  PETSC_LIBDIR
dnl     PETSC_LIB_CONTRIB	  PETSC_LIB_TYPE
dnl     PETSC_LIB_DM		  PETSC_LIB_STATUS
dnl     PETSC_LIB_KSP           Complete Substitutions
dnl     PETSC_LIB_MAT		  PETSC_INC
dnl     PETSC_LIB_SNES		  PETSC_LIBS
dnl     PETSC_LIB_TS
dnl     PETSC_LIB_VEC
dnl
dnl Version:    $Id: tx_petsc3.m4 3773 2010-12-23 00:30:38Z cary $
dnl
dnl Copyright 2001-2010, Tech-X Corporation.  Redistribution allowed
dnl provided this copyright statement remains intact.
dnl
dnl ######################################################################

striplibs() {
  for i in $*; do
    case $i in
      -l*)
        shift
        ;;
      *)
        break
        ;;
    esac
  done
  echo $*
}

dnl Allow for the disabling of petsc
AC_ARG_ENABLE(
  [petsc],
  [AC_HELP_STRING([--disable-petsc],
  [Disables the petsc package from being searched])],
  [tx_petsc_enable=${enableval}])

dnl Look for PETSc directories in the usual suspect locations.  First
dnl check for configure flags, then env variable, then defaults.

if test "x$tx_petsc_enable" != "xno"; then
  dnl AC_MSG_CHECKING([for petsc directory])

  AC_ARG_WITH(
    [petsc-dir],
    [AC_HELP_STRING([--with-petsc-dir=<dir>],
      [Root PETSc directory in which to search])],
    [PETSC_DIR=${withval}])

  if test -n "$PETSC_DIR"; then
    tx_petsc_path=${PETSC_DIR}
  else
    unset tx_petsc_path
    PETSC_SP=$SUPRA_SEARCH_PATH
    if test "$parallel" = yes; then
      for i in `echo $PETSC_SP | tr ':' ' '`; do
        tx_petsc_path="$tx_petsc_path:$i/petsc-par:$i/petscmpi"
      done
    else
      for i in `echo $PETSC_SP | tr ':' ' '`; do
        tx_petsc_path="$tx_petsc_path:$i/petsc"
      done
    fi
  fi

  AC_ARG_WITH(
    [petsc-arch],
    [AC_HELP_STRING(
      [--with-petsc-arch=<arch>],
      [Architecture of the PETSc build])],
    [PETSC_ARCH=${withval}])

  if test "$parallel" = yes; then
    TX_LOCATE_PKG(
      [PETSC],
      [$tx_petsc_path],
      [petsc.h],
      [petscts, petscsnes, petscksp, petscdm, petscmat, petscvec, petsc, dmumps, cmumps, mumps_common, smumps, zmumps, scalapack, blacs, pord, superlu_dist_2.4, HYPRE, parmetis, metis],
      [include],
      [lib])
  else
    TX_LOCATE_PKG(
      [PETSC],
      [$tx_petsc_path],
      [mpif.h,petsc.h],
      [petscts, petscsnes, petscksp, petscdm, petscmat, petscvec, petsc, superlu_4.0],
      [include],
      [lib])
# Add in mpiuni if present.  Do this way or the incdir is set
# to $PETSC_INCDIR/mpiuni and PETSC_DIR is set to the include dir.
    if test -d $PETSC_INCDIR/mpiuni; then
      PETSC_INC="$PETSC_INC -I$PETSC_INCDIR/mpiuni"
    fi
  fi

  echo PETSC_INCDIR = $PETSC_INCDIR
  if test -z "$PETSC_LIBDIR"; then
    dnl AC_MSG_RESULT([not found])
    if test -z "$PETSC_REQUIRED"; then
      AC_MSG_WARN([PETSC directory could not be found!  Use --with-petsc-dir to specify one manually])
    else
      AC_MSG_ERROR([PETSC directory could not be found!  Use --with-petsc-dir to specify one manually])
    fi
  else
    PETSC_DIR=`(cd $PETSC_INCDIR/..; pwd -P)`
    dnl AC_MSG_RESULT([$PETSC_DIR])
  fi
  echo PETSC_DIR = $PETSC_DIR

#
# Get just names and dirs as needed by uedge
#

  if test -n "$PETSC_LIBS"; then
# Remove duplicates and get just incdirs
    unset PETSC_INCDIRS
    petscinc="$PETSC_INC"
    for i in $petscinc; do
      case $i in
        -I*)
          pcinc=`echo $i | sed 's/^-I//'`
          if ! echo $pcinc | egrep -q -- "(^| )$pcinc($| )"; then
            PETSC_INCDIRS="$PETSC_INCDIRS $pcinc"
            PETSC_INC="$PETSC_INC -I$pcinc"
          fi
          ;;
      esac
    done
    TX_PRINT_VAR(PETSC_INCDIRS)
    unset PETSC_LIBDIRS
    unset PETSC_LIBNAMES
    for i in $PETSC_LIBS; do
      case $i in
        -l*)
          pclib=`echo $i | sed 's/^-l//'`
          PETSC_LIBNAMES="$PETSC_LIBNAMES $pclib"
          ;;
        -L*)
          pclibdir=`echo $i | sed 's/^-L//'`
          PETSC_LIBDIRS="$PETSC_LIBDIRS $pclibdir"
          ;;
      esac
    done
    TX_PRINT_VAR(PETSC_LIBDIRS)
    TX_PRINT_VAR(PETSC_LIBNAMES)
  fi
  AC_SUBST(PETSC_INCDIRS)
  AC_SUBST(PETSC_LIBDIRS)
  AC_SUBST(PETSC_LIBNAMES)

#
# Parallel lbraries
# Create petsc makefile for getting mpi includes and directories
#

cat <<MFILE > petscmake$$
PETSC_DIR  = ${PETSC_DIR}
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
MFILE

# Get mpi includes
  tmp=`make -f petscmake$$ getmpiincludedirs`
# Remove mpiuni directory as points to build
  unset petscmpiinc
  for i in $tmp; do
    case $i in
      -I.*mpiuni)
        ;;
      *)
# Add unique include dirs
        if ! echo $petscmpiinc | egrep -q "(^| )$i($| )"; then
          petscmpiinc="$petscmpiinc $i"
        fi
        ;;
    esac
  done
  # echo petscmpiinc = $petscmpiinc
  petscmpiincdirs=`echo $petscmpiinc | sed -e 's/^-I//' -e 's/ -I/ /g'`
  # echo petscmpiincdirs = $petscmpiincdirs
# Some MPI implementations put mpi.mod in the lib directory,
# and this is not reported by petsc.
  havempimod=false
  for j in $petscmpiincdirs; do
    if test -f $j/mpi.mod; then
      havempimod=true
      break
    # else
      # echo mpi.mod not found in $j.
    fi
  done
  # echo havempimod = $havempimod
  if ! $havempimod; then
    echo Looking for mpi.mod.
    for j in $petscmpiincdirs; do
      if test -f $j/../lib/mpi.mod; then
        petscmpimoddir=`(cd $j/../lib; pwd -P)`
        petscmpiinc="$petscmpiinc -I$petscmpimoddir"
        break
      fi
    done
  fi
  PETSC_MPI_INC="$petscmpiinc"
  PETSC_MPI_INCDIRS=`echo $PETSC_MPI_INC | sed -e 's/^-I//g' -e 's/ -I/ /g'`

# Get mpi libraries
  petscalllibs=`make -f petscmake$$ getlinklibs`
  # echo petscalllibs = $petscalllibs
  rm petscmake$$
  # petscremlibs=`echo $petscalllibs | sed "s?^.*$PETSC_LIBDIR??"`
  # petscremlibs=`striplibs $petscremlibs`
# Squeeze out spaces between -L and dir
  # petscremlibs=`echo $petscalllibs | sed -e 's/^.*-lsuperlu_[34].[0-9]//' -e 's/-L  */-L/'`
  petscremlibs=`echo $petscalllibs | sed -e 's/-L  */-L/'`
# Keep one instance of petsclibdir
  petscremlibs="$petscremlibs -L$PETSC_LIBDIR"
  # echo petscremlibs = $petscremlibs

# Correct for framework stuff
  unset petscmpilibs
  unset petscmpilibflags
  framework=false
  for i in $petscremlibs; do
    if $framework; then
      var=`echo $i | sed 's/^-l//'`
      petscmpilibs="$petscmpilibs $var"
      framework=false
      continue
    fi
    unset lib
    unset flag
    case $i in
      -framework)
        petscmpilibs="$petscmpilibs -framework"
        framework=true
        ;;
      -lAccelerate) # Skip this, as can come from only grabnext
        ;;
      -L/usr/lib/gcc/* | -L/usr/libexec/gcc/* | -L/usr/local/libexec/gcc/i386-apple-darwin* | -L/usr/local/lib/gcc/i386-apple-darwin* | -Wl,-rpath,* | -llapack | -lblas |  -lgfortran* | -lstdc++ | -lgcc* | -lSystem | -ldl | -lm)
# Skip system and compiler libs
        ;;
      *hdf* | *netcdf* )
        ;;
      -lpetscts | -lpetscsnes | -lpetscksp | -lpetscdm | -lpetscmat | -lpetscvec | -lpetsc | -lcmumps | -ldmumps | -lsmumps | -lzmumps | -lmumps_common | -lpord | -lscalapack | -lblacs | -lsuperlu_dist_2.4 | -lsuperlu*_[2-4].[0-9] | -lHYPRE | -llapack | -lblas | -lparmetis | -lmetis)
# Skip petsc libs
        ;;
      -L*)
        flag=$i
        ;;
      -lmpif90 | -lmpi_f77 | -lmpi_cxx | -lfmpich | -lmpichcxx | -lmpichf90)
# These go before
        lib=$i
        ;;
      *)
# The rest go after
        lib=$i
        ;;
    esac
# Grab flags
    if test -n "$flag"; then
      if ! echo $petscmpilibflags | egrep -q "(^| )$flag($| )"; then
        petscmpilibflags="$petscmpilibflags $flag"
      fi
    fi
# Grab libs
    if test -n "$lib"; then
      if ! echo $petscmpilibs | egrep -q "(^| )$lib($| )"; then
        petscmpilibs="$petscmpilibs $lib"
      fi
    fi
  done

# Removing duplicates causes mpi_f77 not to precede the other parallel libs,
# so fix
  if echo $petscmpilibs | egrep -q "(^| )-lmpi_f77($| )"; then
    petscmpilibs="-lmpi_f77 "`echo $petscmpilibs | sed -e 's/ -lmpi_f77 / /' -e 's/ -lmpi_f77$//'`
  fi
  petscmpilibs="$petscmpilibflags $petscmpilibs"
  # echo petscmpilibs = $petscmpilibs
  PETSC_MPI_LIBS="$petscmpilibs"
  unset PETSC_MPI_LIBDIRS
  unset PETSC_MPI_LIBNAMES
  for i in $PETSC_MPI_LIBS; do
    case $i in
     -l*)
        libname=`echo $i | sed 's/^-l//'`
        PETSC_MPI_LIBNAMES="$PETSC_MPI_LIBNAMES $libname"
        ;;
     -L*)
        dirname=`echo $i | sed 's/^-L//'`
        PETSC_MPI_LIBDIRS="$PETSC_MPI_LIBDIRS $dirname"
        ;;
    esac
  done

fi

# Reverse default definition, as part is good here
AM_CONDITIONAL(HAVE_PETSC, test -n "$PETSC_LIBS")
if test -n "$PETSC_LIBS"; then
  AC_DEFINE(HAVE_PETSC, [], [Defined if PETSC found])
fi

AC_SUBST(PETSC_DIR)
AC_SUBST(PETSC_ARCH)
TX_PRINT_VAR(PETSC_ARCH)

AC_SUBST(PETSC_MPI_INC)
AC_SUBST(PETSC_MPI_INCDIRS)
TX_CLEAN_LIBS([PETSC_MPI_LIBS])
AC_SUBST(PETSC_MPI_LIBS)
AC_SUBST(PETSC_MPI_RPLIBS)
AC_SUBST(PETSC_MPI_LTLIBS)
AC_SUBST(PETSC_MPI_ALIBS)
AC_SUBST(PETSC_MPI_LIBDIRS)
AC_SUBST(PETSC_MPI_LIBNAMES)

if test "x$tx_petsc_enable" != "xno"; then
  echo >> $config_summary_file
  echo PETSC MPI information found with >> $config_summary_file
  TX_PRINT_VAR(PETSC_MPI_INC)
  TX_PRINT_VAR(PETSC_MPI_INCDIRS)
  TX_PRINT_VAR(PETSC_MPI_LIBS)
  TX_PRINT_VAR(PETSC_MPI_RPLIBS)
  TX_PRINT_VAR(PETSC_MPI_LTLIBS)
  TX_PRINT_VAR(PETSC_MPI_ALIBS)
  TX_PRINT_VAR(PETSC_MPI_LIBDIRS)
  TX_PRINT_VAR(PETSC_MPI_LIBNAMES)
else
  echo >> $config_summary_file
  echo "NOT using petsc" >> $config_summary_file
fi

