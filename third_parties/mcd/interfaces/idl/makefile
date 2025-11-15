#
# Makefile for IDL and MCD
# make -f makefile
#
# Path to (top directory of) NetCDF distribution 
NETCDF_DIR = /path/to/netcdf
# Path to the MCD Fortan routines (MCD.F90)
FROUTINES_DIR = /path/to/MCD.F90
# Fortran compiler (should be the same that was used to compile the NetCDF library)
FCOMPILER = gfortran

call_mcd.so : MCD.o julian.o julian_wrapf.o call_mcd_wrapf.o call_mcd4idl.o
	$(FCOMPILER) -shared -fPIC -o call_mcd.so MCD.o julian.o julian_wrapf.o call_mcd_wrapf.o call_mcd4idl.o $(NETCDF_DIR)/lib/libnetcdff.so

MCD.o : $(FROUTINES_DIR)/MCD.F90
	$(FCOMPILER) -shared -fPIC -c $(FROUTINES_DIR)/MCD.F90 -I$(NETCDF_DIR)/include
call_mcd_wrapf.o : call_mcd_wrapf.f90
	$(FCOMPILER) -shared -fPIC -c call_mcd_wrapf.f90
call_mcd4idl.o   : call_mcd4idl.f90
	$(FCOMPILER) -shared -fPIC -c call_mcd4idl.f90
julian.o : julian.F90        
	$(FCOMPILER) -shared -fPIC -c julian.F90
julian_wrapf.o   : julian.o julian_wrapf.f90 
	$(FCOMPILER) -shared -fPIC -c julian_wrapf.f90

