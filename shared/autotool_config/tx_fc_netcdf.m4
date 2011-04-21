#
# SYNOPSIS
#
#   TX_FC_ALLOCATABLE_COMPONENT
#
# DESCRIPTION
#
#   Check that netcdf.mod is compatible with fortran compiler
#
#   $Id: tx_fc_netcdf.m4 3507 2010-05-04 12:32:14Z cary $
#
AC_DEFUN([TX_FC_NETCDF],[
AC_CACHE_CHECK([whether netcdf module works],
ac_cv_netcdf_mod_is_ok,
[AC_LANG_PUSH(Fortran)
i=0
while test \( -f tmpdir_$i \) -o \( -d tmpdir_$i \) ; do
  i=`expr $i + 1`
  echo $i
done
mkdir tmpdir_$i
cd tmpdir_$i
AC_COMPILE_IFELSE([
!234567
	program test
	use netcdf
	implicit none
	integer :: ier, ncid
	ier = nf90_open('/tmp/foo.nc', NF90_NOWRITE, ncid)
	end program test
],
[ac_cv_netcdf_mod_is_ok="yes"],
[ac_cv_netcdf_mod_is_ok="no"])
cd ..
rm -fr tmpdir_$i
AC_LANG_POP(Fortran)
])
if test "$ac_cv_netcdf_mod_is_ok" = yes; then
   AC_DEFINE(HAVE_FC_NETCDF_MOD,,[define if netcdf module works])
fi
TX_FC_NETCDF_MOD_IS_OK="$ac_cv_netcdf_mod_is_ok"
AC_SUBST(TX_FC_NETCDF_MOD_IS_OK)
])

