# \file CMakeLists.txt
 #    This is the CMake file for the Input sub-directory.
 #
 #    Path              : /Astrodynamics/Input/
 #    Version           : 3
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
 #      111026    K. Kumar          Adapted file so all headers show in project tree in Qt Creator.
 #/

# Add source files.
set(INPUT_SOURCES
  "${SRCROOT}${INPUTDIR}/fileReader.cpp"
  "${SRCROOT}${INPUTDIR}/textFileReader.cpp"
  "${SRCROOT}${INPUTDIR}/twoLineElementData.cpp"
  "${SRCROOT}${INPUTDIR}/twoLineElementsTextFileReader.cpp"
)

# Add header files.
set(INPUT_HEADERS 
  "${SRCROOT}${INPUTDIR}/fileReader.h"
  "${SRCROOT}${INPUTDIR}/textFileReader.h"
  "${SRCROOT}${INPUTDIR}/twoLineElementData.h"
  "${SRCROOT}${INPUTDIR}/twoLineElementsTextFileReader.h"
)

# Add static libraries.
add_library(tudat_input STATIC ${INPUT_SOURCES} ${INPUT_HEADERS})
setup_tudat_library_target(tudat_input "${SRCROOT}${INPUTDIR}")
setup_headers("${INPUT_HEADERS}" "${INPUTDIR}")

# Add unit tests.
add_executable(test_TextFileReader "${SRCROOT}${INPUTDIR}/unitTestTextFileReader.cpp")
setup_custom_test_program(test_TextFileReader "${SRCROOT}${INPUTDIR}")
target_link_libraries(test_TextFileReader tudat_astrodynamics tudat_input tudat_output tudat_mathematics tudat_basics ${Boost_LIBRARIES})

add_executable(test_TwoLineElementsTextFileReader "${SRCROOT}${INPUTDIR}/unitTestTwoLineElementsTextFileReader.cpp")
setup_custom_test_program(test_TwoLineElementsTextFileReader "${SRCROOT}${INPUTDIR}")
target_link_libraries(test_TwoLineElementsTextFileReader tudat_input tudat_output tudat_astrodynamics tudat_mathematics tudat_basics ${Boost_LIBRARIES})

# End of file.
