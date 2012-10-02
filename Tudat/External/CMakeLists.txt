 #    Copyright (c) 2010-2012 Delft University of Technology.	
 #	
 #    This software is protected by national and international copyright.	
 #    Any unauthorized use, reproduction or modification is unlawful and	
 #    will be prosecuted. Commercial and non-private application of the	
 #    software in any form is strictly prohibited unless otherwise granted	
 #    by the authors.	
 #	
 #    The code is provided without any warranty; without even the implied	
 #    warranty of merchantibility or fitness for a particular purpose.	
 #	
 #    Changelog	
 #      YYMMDD    Author            Comment	
 #      110820    S.M. Persson      File created.	
 #      111025    K. Kumar          Adapted file to work with Revision 194.	
 #      111026    K. Kumar          Adapted file so all headers show in project tree in Qt Creator.	
 #      120202    K. Kumar          Removed astrodynamicsFunctions.h/.cpp and unit test	
 #      120717    D. Dirkx          Modified for SpiceInterface.	
 #		

 # Define the main sub-directories.

if(USE_CSPICE)
  set(SPICEINTERFACEDIR "${EXTERNALDIR}/SpiceInterface")

  # Add subdirectories.
  add_subdirectory("${SRCROOT}${SPICEINTERFACEDIR}")

  # Get target properties for static libraries.
  get_target_property(SPICEINTERFACESOURCES tudat_spice_interface SOURCES)

  # Set SPICE interface directory.
  set(SPICEINTERFACEDIR "${EXTERNALDIR}/SpiceInterface")
endif( )
