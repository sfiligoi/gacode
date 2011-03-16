dnl ######################################################################
dnl
dnl File:	fftw.m4
dnl
dnl Purpose:	Looks for single and double precision versions of 
dnl		FFTW libraries
dnl
dnl Version:	$Id: fftwps.m4 3464 2010-04-05 12:42:26Z cary $
dnl
dnl Tech-X configure system
dnl
dnl ######################################################################

dnl ######################################################################
dnl
dnl Allow the user to specify an overall fftw directory.  If specified,
dnl we look for include and lib under this.
dnl
dnl ######################################################################

AC_ARG_WITH(fftw-dir,[  --with-fftw-dir=<location of fftw installation> ],FFTW_DIR="$withval",FFTW_DIR="")
