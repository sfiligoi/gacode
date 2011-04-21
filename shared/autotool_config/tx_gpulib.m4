dnl ######################################################################
dnl
dnl File:	tx_gpulib.m4
dnl
dnl Purpose:	Determine the locations of gpulib includes and libraries
dnl
dnl Version: $Id: tx_gpulib.m4 3743 2010-11-25 16:25:37Z cary $
dnl
dnl Tech-X configure system
dnl
dnl ######################################################################

dnl Includes functions from txsearch.m4
builtin(include, config/txsearch.m4)

dnl We will enable GPLib by default. Use --disable-gpulib to disable.

AC_ARG_ENABLE( [gpulib], [AC_HELP_STRING([--disable-gpulib],
	[Disables the inclusion of GPU-based computation])],
	[tx_gpulib_enable=${enableval}], [tx_gpulib_enable=yes])

AC_ARG_ENABLE( [gpulib-emulator], [AC_HELP_STRING([--enable-gpulib-emulator],
	[Enables the use of gpu-emulator-based computation])],
	[tx_gpulib_enable_emulator=${enableval}], [tx_gpulib_enable_emulator=no])

# If emulator is enabled, then disable gpulib.

if test "$tx_gplulib_enable_emulator" = yes; then
  tx_gpulib_enable=no
fi

AC_ARG_WITH(gpulib-dir, AC_HELP_STRING([--with-gpulib-dir],
        [gpulib installation directory]),
        [GPULIB_DIR="$withval"; GPULIB_PATH=$GPULIB_DIR])

if test -n "$GPULIB_DIR"; then
  GPULIB_SP=$GPULIB_DIR
else
  GPULIB_SP=$SUPRA_SEARCH_PATH
fi

for dir in `echo $GPULIB_SP | tr ':' ' '`; do
  GPULIB_PATH="$GPULIB_PATH:$dir/gpulib"
done

AC_ARG_WITH(cuda-dir, AC_HELP_STRING([--with-cuda-dir],
        [cuda installation directory]),
        [CUDA_DIR="$withval"; CUDA_PATH=$CUDA_DIR])


# HACK: CUDA isn't likely to be in SUPRA_SEARCH_PATH, so we add /usr/local

if test -n "$CUDA_DIR"; then
  CUDA_SP=$CUDA_DIR
else
  CUDA_SP=/usr/local:$SUPRA_SEARCH_PATH
fi

for dir in `echo $CUDA_SP | tr ':' ' '`; do
  CUDA_PATH="$CUDA_PATH:$dir/cuda"
done


AC_ARG_WITH(thrust-dir, AC_HELP_STRING([--with-thrust-dir],
	[thrust installation directory]),
	[THRUST_DIR="$withval"; THRUST_PATH=$THRUST_DIR])

if test -n "$THRUST_DIR"; then
  THRUST_PATH=$THRUST_DIR
else
  THRUST_PATH=$SUPRA_SEARCH_PATH
fi


dnl ######################################################################
dnl
dnl Find
dnl
dnl ######################################################################

if test $tx_gpulib_enable = yes; then
  TX_LOCATE_PKG(
	[GPULIB],
	[$GPULIB_PATH],
	[gpuVectorOp.h, gpuPhysicsOp.h],
	[gpulib],
	[vectorOp, physicsOp],
	[])


  if test -n "$GPULIB_LIBS"; then
    TX_LOCATE_PKG(
	[CUDA],
	[$CUDA_PATH],
	[cuda.h, cublas.h, cuda_runtime.h, driver_types.h, cuda_runtime_api.h],
	[cudart, cublas],
	[include],
	[lib64])
  fi

  if test -n "$GPULIB_LIBS"; then
    TX_LOCATE_PKG(
	[THRUST],
	[$THRUST_PATH, $GPULIB_PATH],
	[device_vector.h],
	[-],
	[thrust],
	[-])
  fi

fi

if test $tx_gpulib_enable_emulator = yes; then
   TX_LOCATE_PKG(
	[GPULIB],
	[$GPULIB_PATH],
	[gpuVectorOp.h,gpuPhysicsOp.h],
	[GPULibemu],
	[include],
	[lib])

   TX_LOCATE_PKG(
	[CUDA],
	[$CUDA_PATH],
	[cuda.h,cublas.h,cuda_runtime.h,driver_types.h],
	[cudart,cublasemu],
	[include],
	[lib])

fi

if test -n "$GPULIB_LIBS"; then
  ac_cv_have_gpulib=yes
else
  ac_cv_have_gpulib=no
fi

AM_CONDITIONAL(HAVE_GPULIB, test -n "$GPULIB_LIBS")
AM_CONDITIONAL(HAVE_THRUST, test -n "$THRUST_INC")

