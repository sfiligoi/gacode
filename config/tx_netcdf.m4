dnl ######################################################################
dnl
dnl File:	tx_netcdf.m4
dnl
dnl Purpose:	Find the various libraries that are installed by the
dnl               netcdf package.
dnl
dnl Version: $Id: tx_netcdf.m4 3808 2011-03-09 03:17:54Z kruger $
dnl
dnl Tech-X configure system.  Copyright 2007-2010.  Freely redistributable
dnl provided this copyright remains intact.
dnl
dnl ######################################################################

NETCDF_SP=$SUPRA_SEARCH_PATH
unset NETCDF_PATH
# echo back_end_node = $back_end_node
# echo parallel = $parallel
if test -n "$back_end_node" -a "$parallel" = yes; then
  for i in `echo $NETCDF_SP | tr ':' ' '`; do
    NETCDF_PATH="$NETCDF_PATH:$i/netcdf-ben"
  done
fi
if test "$parallel" = yes; then
  for i in `echo $NETCDF_SP | tr ':' ' '`; do
    NETCDF_PATH="$NETCDF_PATH:$i/netcdf-par"
  done
fi
for i in `echo $NETCDF_SP | tr ':' ' '`; do
  NETCDF_PATH="$NETCDF_PATH:$i/netcdf"
done
for i in `echo $NETCDF_SP | tr ':' ' '`; do
  NETCDF_PATH="$NETCDF_PATH:$i"
done

case "$TX_FORTRAN_MODCAP" in
  ucname-*)
    TX_LOCATE_PKG(
      [NETCDF],
      [$NETCDF_PATH],
      [netcdf.h, NETCDF.mod],
      [netcdff, netcdf])
    ;;
  lcname-*)
    TX_LOCATE_PKG(
      [NETCDF],
      [$NETCDF_PATH],
      [netcdf.h, netcdf.mod],
      [netcdff, netcdf])
    ;;
esac

# Find any libraries is good
AM_CONDITIONAL(HAVE_NETCDF, test -n "$NETCDF_LIBS")

