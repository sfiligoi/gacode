#---------------------------------------------------
# JAC {Portland Group}  
#---------------------------------------------------

# Compilers and flags

FC	= mpif90 -module $(GACODE_ROOT)/modules
F77	= mpif90

FMATH	= -r8 
FOPT	= -fastsse 
FDEBUG	= -g -Ktrap=fp -Mbounds -Minform=warn
FBOUND	= -Minform=warn -Mbounds
FPROF	= -Mprof=lines

# System math libraries

#LMATH = /usr/lib/liblapack.a /usr/lib/libblas.a
#LMATH = /usr/local/depot/PGI/pgi-7.0/linux86/7.0-7/lib/liblapack.a \
#        /usr/local/depot/PGI/pgi-7.0/linux86/7.0-7/lib/libblas.a
LMATH = /usr/local/depot/PGI/pgi-12.4/linux86/12.4/lib/liblapack.a \
        /usr/local/depot/PGI/pgi-12.4/linux86/12.4/lib/libblas.a

# netCDF (not optional any more)

NETCDF = /usr/local/depot/netcdf-4.0.1/lib/libnetcdf.a
NETCDF_INC = /usr/local/depot/netcdf-4.0.1/include

# Archive 

ARCH = ar cr

# The following line _might_ be needed if you get the error message about m2c being unable to compile the .mod files. However, it is more likely that you mod files are not stored at the correct location. 
#%.o : %.mod
