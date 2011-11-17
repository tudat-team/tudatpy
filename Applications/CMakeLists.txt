# \file CMakeLists.txt
 #    This is the CMake file for the Applications sub-directory.
 #
 #    Path              : /Applications/
 #    Version           : 2
 #    Check status      : Checked
 #
 #    Author            : S.M. Persson
 #    Affiliation       : McGill University, Department of Mechanical Engineering
 #    E-mail address    : mikael.s.persson@gmail.com
 #
 #    Author            : K. Kumar
 #    Affiliation       : Delft University of Technology
 #    E-mail address    : K.Kumar@tudelft.nl
 #
 #    Author            : S.M. Persson
 #    Affiliation       : McGill University, Department of Mechanical Engineering
 #    E-mail address    : mikael.s.persson@gmail.com
 #
 #    Date created      : 20 August, 2011
 #    Last modified     : 25 October, 2011
 #
 #    References
 #
 #    Notes
 #
 #    Copyright (c) 2010-2011 Delft University of Technology.
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
 #/

# Set source files.
#set(APPLICATIONS_SOURCES 
#  ""
#)

# Set header files.
#set(APPLICATIONS_HEADERS 
#  ""
#)

# Add as static library.
#add_library(tudat_applications STATIC ${APPLICATIONS_SOURCES})
#setup_tudat_library_target(tudat_applications "${SRCROOT}${APPLICATIONSDIR}")
#setup_headers("${APPLICATIONS_HEADERS}" "${APPLICATIONSDIR}")

add_executable(application_ExampleEarthOrbitingSatellite "${SRCROOT}${APPLICATIONSDIR}/exampleEarthOrbitingSatellite.cpp")
setup_tudat_executable_target(application_ExampleEarthOrbitingSatellite "${SRCROOT}${APPLICATIONSDIR}")
target_link_libraries(application_ExampleEarthOrbitingSatellite tudat_astrodynamics tudat_input tudat_output tudat_mathematics tudat_basics ${Boost_LIBRARIES})

# End of file.
