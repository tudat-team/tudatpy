# For MATLAB access to MCD.F90
# MCD.o 
#
# Path to (top directory of) NetCDF distribution 
NETCDF_DIR = /path/to/netcdf
# Path to the Fortan routines (MCD.F90)
FROUTINES_DIR = ..
# Fortran compiler (should be the same that was used to compiler NetCDF)
FCOMPILER = gfortran

prog: mcd.o

mcd.o: $(FROUTINES_DIR)/MCD.F90
	$(FCOMPILER) -c -fPIC -I$(NETCDF_DIR)/include $(FROUTINES_DIR)/MCD.F90
