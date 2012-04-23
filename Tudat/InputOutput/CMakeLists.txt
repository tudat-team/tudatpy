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
 #      120202    B. Tong Minh      File copied from Tudat core
 #

# Add source files.
set(INPUTOUTPUT_SOURCES
  "${SRCROOT}${INPUTOUTPUTDIR}/basicInputOutput.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/twoLineElementData.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/twoLineElementsTextFileReader.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/matrixTextFileReader.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/fieldValue.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/parsedDataVectorUtilities.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/textParser.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/whiteSpaceParser.cpp"
)

# Add header files.
set(INPUTOUTPUT_HEADERS 
  "${SRCROOT}${INPUTOUTPUTDIR}/basicInputOutput.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/twoLineElementData.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/twoLineElementsTextFileReader.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/matrixTextFileReader.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/fieldTransform.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/fieldType.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/fieldValue.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/parsedDataVectorUtilities.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/parser.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/textParser.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/whiteSpaceParser.h"
)

# Add static libraries.
add_library(tudat_input_output STATIC ${INPUTOUTPUT_SOURCES} ${INPUTOUTPUT_HEADERS})
setup_tudat_library_target(tudat_input_output "${SRCROOT}${INPUTOUTPUTDIR}")
setup_headers("${INPUTOUTPUT_HEADERS}" "${INPUTOUTPUTDIR}")

# Add unit tests.
add_executable(test_MatrixTextFileReader "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestMatrixTextFileReader.cpp")
setup_custom_test_program(test_MatrixTextFileReader "${SRCROOT}${INPUTDIR}")
target_link_libraries(test_MatrixTextFileReader tudat_input_output ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_TwoLineElementsTextFileReader "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestTwoLineElementsTextFileReader.cpp")
setup_custom_test_program(test_TwoLineElementsTextFileReader "${SRCROOT}${INPUTDIR}")
target_link_libraries(test_TwoLineElementsTextFileReader tudat_input_output ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_BasicInputOutput "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestBasicInputOutput.cpp")
setup_custom_test_program(test_BasicInputOutput "${SRCROOT}${INPUTDIR}")
target_link_libraries(test_BasicInputOutput tudat_input_output ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_WhiteSpaceParser "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestWhitespaceParser.cpp")
setup_custom_test_program(test_WhiteSpaceParser "${SRCROOT}${INPUTDIR}")
target_link_libraries(test_WhiteSpaceParser tudat_input_output ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})
