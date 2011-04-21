#
# SYNOPSIS
#
#   TX_FC_HAS_INTRINSIC
#
# DESCRIPTION
#
#  Check whether compiler supports supplied subroutine as intrinsic
#
#  $Id: tx_fc_has_intrinsic.m4 3464 2010-04-05 12:42:26Z cary $
#
AC_DEFUN([TX_FC_HAS_INTRINSIC],[
AS_VAR_PUSHDEF([type_var], [ac_cv_has_intrinsic_$1])
AC_CACHE_CHECK([whether compiler has intrinsic subroutine $1],
ac_cv_has_intrinsic_$1,
[AC_LANG_PUSH(Fortran)
i=0
while test \( -f tmpdir_$i \) -o \( -d tmpdir_$i \) ; do
  i=`expr $i + 1`
done
mkdir tmpdir_$i
cd tmpdir_$i
AC_LINK_IFELSE([
!234567
      program test
       call $1
      end
  ],
  [ac_cv_has_intrinsic_$1="yes"],
  [ac_cv_has_intrinsic_$1="no"])
cd ..
rm -fr tmpdir_$i
AC_LANG_POP(Fortran)
])
if test "$ac_cv_has_intrinsic_$1" = yes; then
   AC_DEFINE(HAVE_FC_INTRINSIC_$1,,[define if fortran compiler has intrinsic subroutine $1])
fi
TX_FORTRAN_HAS_INTRINSIC_$1="$ac_cv_has_intrinsic_$1"

dnl ######################################################################
dnl
dnl Print out configuration information
dnl
dnl ######################################################################

if test -n "$config_summary_file"; then
   echo                                      >> $config_summary_file
   if test "$TX_FORTRAN_HAS_INTRINSIC_$1" = yes; then
      echo "Fortran compiler supports intrinsic $1"  >> $config_summary_file
   else
      echo "Fortran compiler does NOT support instrinsic $1" >> $config_summary_file
   fi
fi
AC_SUBST(TX_FORTRAN_HAS_INTRINSIC_$1)
])

