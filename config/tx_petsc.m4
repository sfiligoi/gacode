dnl ######################################################################
dnl
dnl File:     tx_petsc2.m4
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
dnl Version:    $Id: tx_petsc.m4 3620 2010-08-18 14:59:25Z cary $
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
    tx_petsc_path=$PETSC_DIR
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
      [petscts, petscsnes, petscksp, petscdm, petscmat, petscvec, petsc, cmumps, dmumps, smumps, zmumps, mumps_common, pord, scalapack, blacs, superlu_dist_2.3, HYPRE, parmetis, metis],
      [include],
      [lib])
  else
    TX_LOCATE_PKG(
      [PETSC],
      [$tx_petsc_path],
      [mpif.h, petsc.h],
      [petscts, petscsnes, petscksp, petscdm, petscmat, petscvec, petsc, mpiuni, superlu_3.1],
      [include, include/mpiuni],
      [lib])
  fi

  echo PETSC_INCDIR = $PETSC_INCDIR
  if test -z "$PETSC_LIBDIR"; then
    dnl AC_MSG_RESULT([not found])
    AC_MSG_ERROR([PETSC directory could not be found!  Use --with-petsc-dir to specify one manually])
  else
    PETSC_DIR=`(cd $PETSC_INCDIR/..; pwd -P)`
    dnl AC_MSG_RESULT([$PETSC_DIR])
  fi
  echo PETSC_DIR = $PETSC_DIR

# Create petsc makefile for getting mpi includes and directories
  if test -f "${PETSC_DIR}/conf/base"; then
    tx_petsc_conf="${PETSC_DIR}/conf/base"
  elif test -f "${PETSC_DIR}/bmake/common/base"; then
    tx_petsc_conf="${PETSC_DIR}/bmake/common/base"
  else
    AC_MSG_ERROR([petsc configure summary was not found or is not a regular file.  Use --with-petsc-dir to specify a directory])
  fi
cat <<MFILE > petscmake$$
PETSC_DIR  = ${PETSC_DIR}
include $tx_petsc_conf
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
        petscmpiinc="$petscmpiinc $i"
        ;;
    esac
  done
  # echo petscmpiinc = $petscmpiinc
  petscmpiincdir=`echo $petscmpiinc | sed -e 's/^-I//'`
  # echo petscmpiincdir = $petscmpiincdir
# Some MPI implementations put mpi.mod in the lib directory
  if test -f $petscmpiincdir/../lib/mpi.mod; then
    petscmpimoddir=`(cd $petscmpiincdir/../lib; pwd -P)`
    petscmpiinc="$petscmpiinc -I$petscmpimoddir"
  fi
  PETSC_MPI_INC="$petscmpiinc"

# Get mpi libraries
  petscalllibs=`make -f petscmake$$ getlinklibs`
  # echo petscalllibs = $petscalllibs
  rm petscmake$$
  # petscremlibs=`echo $petscalllibs | sed "s?^.*$PETSC_LIBDIR??"`
  # petscremlibs=`striplibs $petscremlibs`
# Squeeze out spaces between -L and dir
  petscremlibs=`echo $petscalllibs | sed -e 's/^.*-lsuperlu_3.1//' -e 's/-L  */-L/'`
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
      -lpetscts | -lpetscsnes | -lpetscksp | -lpetscdm | -lpetscmat | -lpetscvec | -lpetsc | -lcmumps | -ldmumps | -lsmumps | -lzmumps | -lmumps_common | -lpord | -lscalapack | -lblacs | -lsuperlu_dist_2.3 | -lHYPRE | -llapack | -lblas | -lparmetis | -lmetis)
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

  petscmpilibs="$petscmpilibflags $petscmpilibs"
  # echo petscmpilibs = $petscmpilibs
  PETSC_MPI_LIBS="$petscmpilibs"

fi

AC_SUBST(PETSC_DIR)
AC_SUBST(PETSC_ARCH)
TX_PRINT_VAR(PETSC_ARCH)

AC_SUBST(PETSC_MPI_INC)
TX_CLEAN_LIBS([PETSC_MPI_LIBS])
AC_SUBST(PETSC_MPI_LIBS)
AC_SUBST(PETSC_MPI_RPLIBS)
AC_SUBST(PETSC_MPI_LTLIBS)
AC_SUBST(PETSC_MPI_ALIBS)

if test "x$tx_petsc_enable" != "xno"; then
  echo >> $config_summary_file
  echo PETSC MPI information found with >> $config_summary_file
  TX_PRINT_VAR(PETSC_MPI_INC)
  TX_PRINT_VAR(PETSC_MPI_LIBS)
  TX_PRINT_VAR(PETSC_MPI_RPLIBS)
  TX_PRINT_VAR(PETSC_MPI_LTLIBS)
  TX_PRINT_VAR(PETSC_MPI_ALIBS)
else
  echo >> $config_summary_file
  echo "NOT using petsc" >> $config_summary_file
fi

