#---------------------------------------------------
# MINT 2-core Workstation
#---------------------------------------------------

IDENTITY="2-core MINT workstation"
CORES_PER_NODE=2
NUMAS_PER_NODE=1

FFTWDIR=/usr/lib/x86_64-linux-gnu
ROOT=/home/candy

# Compilers

FC	= mpif90 -std=f2008 -I $(GACODE_ROOT)/modules -I /usr/include
F77	= mpif77
FOMP    = -fopenmp
FMATH	= -fdefault-real-8 -J $(GACODE_ROOT)/modules
FOPT	= -Ofast
FDEBUG	= -Wall -fcheck=all -Wunderflow -finit-real=nan -fbounds-check 

# System math libraries
LMATH = ${ROOT}/OpenBLAS/libopenblas.a ${FFTWDIR}/libfftw3.a ${FFTWDIR}/libfftw3_omp.a

# OPTIONAL NetCDF:
NETCDF=-L/usr/lib -lnetcdff -lnetcdf
NETCDF_INC=/usr/include

# Archive 

ARCH = ar cr
