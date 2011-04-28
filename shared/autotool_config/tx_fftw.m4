dnl ######################################################################
dnl
dnl File:	fftw.m4
dnl
dnl Purpose:	Looks for single and double precision versions of
dnl		FFTW libraries
dnl
dnl Version:	$Id: tx_fftw.m4 3775 2011-01-12 00:33:31Z kruger $
dnl
dnl Tech-X configure system
dnl
dnl ######################################################################

FFTW_SP=$SUPRA_SEARCH_PATH
unset FFTW_PATH
if test -n "$parallel"; then
  for i in `echo $FFTW_SP | tr ':' ' '`; do
    FFTW_PATH="$FFTW_PATH:$i/fftw-par"
  done
fi
for i in `echo $FFTW_SP | tr ':' ' '`; do
  FFTW_PATH="$FFTW_PATH:$i"
done

if test -n "$parallel"; then
    TX_LOCATE_PKG(
      [FFTW],
      [$FFTW_PATH],
      [fftw.h, rfftw.h, fftw_mpi.h, rfftw_mpi.h],
      [fftw,rfftw,fftw_mpi,rfftw_mpi])
else
    TX_LOCATE_PKG(
      [FFTW],
      [$FFTW_PATH],
      [fftw.h, rfftw.h],
      [fftw,rfftw])
fi

# Find any libraries is good
AM_CONDITIONAL(HAVE_FFTW, test -n "$FFTW_LIBS")
