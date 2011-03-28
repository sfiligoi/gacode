dnl###################################################################
dnl
dnl File:       mpilauncher.m4
dnl
dnl Purpose: 	Determine mpilauncher and additional options based on
dnl		platform.
dnl
dnl Version:    $Id: mpilauncher.m4 3656 2010-09-05 23:59:17Z cary $
dnl
dnl###################################################################

smartwhich() {
  count=0
  res=0

  for dir in `echo $PATH | sed 's/:/ /g'`; do
    count=`expr $count + 1`
    if test -x $dir/$1; then
       res=$count
       break
    fi
  done

  echo $res
}

# What launches mpi programs?
AC_MSG_CHECKING(for mpi launcher)
AC_ARG_WITH(mpi-launcher,
    AC_HELP_STRING([--with-mpi-launcher=<mpilauncher>],
        [to specify the binary for launching mpi programs]),
    MPILAUNCHER="$withval")

if test -z "$MPILAUNCHER"; then
# jrc: simpler alternative?
  AC_PATH_PROGS(MPILAUNCHER, aprun mpiexec mpirun, "", $PATH)
else
  MPILAUNCHER=`which $MPILAUNCHER`
fi
AC_MSG_RESULT($MPILAUNCHER)

if test -z "$MPILAUNCHER"; then
  AC_MSG_WARN([Cannot find an mpilauncher (aprun, mpiexec, mpirun) in your path.  Try --with-mpi-launcher=MPILAUNCHER if not on Blue Gene P.])
fi

##########
#
# Add in options for various flavors of mpi
#
##########

MPILAUNCHERBASE=`basename $MPILAUNCHER`
# echo MPILAUNCHERBASE = $MPILAUNCHERBASE
case $MPILAUNCHERBASE in

  aprun)
    ;;

  mpiexec | mpirun)
    # AC_MSG_NOTICE(Using mpiexec or mpirun)
    isopenmpi=`$MPILAUNCHER -V 2>&1 | grep OpenRTE`
    # AC_MSG_NOTICE(isopenmpi = $isopenmpi)
    if test -n "$isopenmpi"; then
      AC_MSG_NOTICE(Using OpenMPI)
      MPIHOSTARGNAME="--hostfile"
      openmpibindir=`dirname $MPILAUNCHER`
      openmpiexec=$openmpibindir/mpiexec
      openmpidir=`dirname $openmpibindir`
      echo openmpidir = $openmpidir
      OMPIVER=`echo $isopenmpi | sed -e 's/^.*) //'`
      # echo "Using openmpi, version = '$OMPIVER'"
# JRC 27feb10: No longer needed with 1.4.1.
      case $OMPIVER in
        1.2.* | 1.3.* | 1.4.0)
          AC_MSG_WARN(OpenMPI version $OMPIVER is broken.  See https://svn.open-mpi.org/trac/ompi/ticket/2043.  You must upgrade.)
          MPIOPTS="-x PYTHONPATH -x SIDL_DLL_PATH -mca btl ^sm"
          defhostfile=`(cd $openmpidir; pwd -P)`/etc/openmpi-default-hostfile
          nontriv=`sed -e '/^ *#/d' -e '/^ *$/d' <$defhostfile`
          if test -n "$nontriv"; then
            MPIOPTS="$MPIOPTS --default-hostfile $defhostfile"
          fi
          ;;
        1.4.[[1-9]] | 1.5.*)
          AC_MSG_NOTICE(OpenMPI version is after 1.4.1)
          MPIOPTS="-x PYTHONPATH -x SIDL_DLL_PATH"
          defhostfile=`(cd $openmpidir; pwd -P)`/etc/openmpi-default-hostfile
          nontriv=`sed -e '/^ *#/d' -e '/^ *$/d' <$defhostfile`
          if test -n "$nontriv"; then
            MPIOPTS="$MPIOPTS --default-hostfile $defhostfile"
          fi
          ;;
        *)
          echo "Version did not match."
          ;;
      esac
    else
      MPIHOSTARGNAME="-machinefile"
      hasenvall=`$MPILAUNCHER --help 2>&1 | grep envall`
      if test -n "$hasenvall"; then
        MPIOPTS="-envall"
      fi
    fi
    ;;

esac
echo MPIOPTS = $MPIOPTS

##########
#
# Keep host file separate, as script may want to choose how to add this info
#
##########

# Construct nodes arg if present
if test -n "$MPIHOSTARGNAME"; then
  if test -n "$PBS_NODEFILE"; then
    MIPHOSTARGS="$MPIHOSTARGNAME $PBS_NODEFILE"
  else
    for i in .. ../.. ../../..; do
      if test -f $i/nodes; then
        nddir=`(cd $i; pwd -P)`
        MIPHOSTARGS="$MPIHOSTARGNAME $nddir/nodes"
        break
      fi
    done
  fi
fi

##########
#
# Print results
#
##########

# echo Printing MPI results.
if test -n "$config_summary_file"; then
  printf "\n" >> $config_summary_file
  printf "Using MPI with:\n" >> $config_summary_file
  TX_PRINT_VAR(MPILAUNCHER)
  TX_PRINT_VAR(MPIOPTS)
  TX_PRINT_VAR(MPIHOSTARGNAME)
  TX_PRINT_VAR(MPIHOSTARGS)
fi

##########
#
# Save for substitutions
#
##########

AC_SUBST(MPILAUNCHER)
AC_SUBST(MPIOPTS)
AC_SUBST(MPIHOSTARGNAME)
AC_SUBST(MPIHOSTARGS)

