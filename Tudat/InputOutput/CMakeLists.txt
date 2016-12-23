 #    Copyright (c) 2010-2015, Delft University of Technology
 #    All rights reserved.
 #
 #    Redistribution and use in source and binary forms, with or without modification, are
 #    permitted provided that the following conditions are met:
 #      - Redistributions of source code must retain the above copyright notice, this list of
 #        conditions and the following disclaimer.
 #      - Redistributions in binary form must reproduce the above copyright notice, this list of
 #        conditions and the following disclaimer in the documentation and/or other materials
 #        provided with the distribution.
 #      - Neither the name of the Delft University of Technology nor the names of its contributors
 #        may be used to endorse or promote products derived from this software without specific
 #        prior written permission.
 #
 #    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 #    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 #    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 #    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 #    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 #    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 #    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 #    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 #    OF THE POSSIBILITY OF SUCH DAMAGE.
 #
 #    Changelog
 #      YYMMDD    Author            Comment
 #      120202    B. Tong Minh      File copied from Tudat core.
 #      130111    K. Kumar          Removed matrix text reader files (ported to Tudat Core).
 #
 #    References
 #   
 #    Notes
 #

# Add source files.
set(INPUTOUTPUT_SOURCES
  "${SRCROOT}${INPUTOUTPUTDIR}/basicInputOutput.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/dictionaryComparer.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/dictionaryTools.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/fieldValue.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/fixedWidthParser.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/linearFieldTransform.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/missileDatcomData.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/missileDatcomReader.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/parsedDataVectorUtilities.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/separatedParser.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/textParser.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/twoLineElementData.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/twoLineElementsTextFileReader.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/matrixTextFileReader.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/streamFilters.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/parseSolarActivityData.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/extractSolarActivityData.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/solarActivityData.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/multiDimensionalArrayReader.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/aerodynamicCoefficientReader.cpp"
)

# Add header files.
set(INPUTOUTPUT_HEADERS 
  "${SRCROOT}${INPUTOUTPUTDIR}/basicInputOutput.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/dictionaryComparer.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/dictionaryEntry.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/dictionaryTools.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/extractor.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/fieldTransform.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/fieldType.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/fieldValue.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/fixedWidthParser.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/linearFieldTransform.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/missileDatcomData.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/missileDatcomReader.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/parsedDataVectorUtilities.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/parser.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/separatedParser.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/textParser.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/twoLineElementData.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/twoLineElementsTextFileReader.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/basicInputOutput.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/matrixTextFileReader.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/streamFilters.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/parseSolarActivityData.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/extractSolarActivityData.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/solarActivityData.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/multiDimensionalArrayReader.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/aerodynamicCoefficientReader.h"
)

# Add unit test files.
set(INPUTOUTPUT_UNITTESTS
  "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestBasicInputOutput.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestMatrixTextFileReader.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestStreamFilters.cpp"
)

# Add static libraries.
add_library(tudat_input_output STATIC ${INPUTOUTPUT_SOURCES} ${INPUTOUTPUT_HEADERS})
setup_tudat_library_target(tudat_input_output "${SRCROOT}${INPUTOUTPUTDIR}")

# Add unit tests.
add_executable(test_MatrixTextFileReader "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestMatrixTextFileReader.cpp")
setup_custom_test_program(test_MatrixTextFileReader "${SRCROOT}${INPUTOUTPUTDIR}")
target_link_libraries(test_MatrixTextFileReader tudat_input_output tudat_basic_astrodynamics ${Boost_LIBRARIES})

add_executable(test_StreamFilters "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestStreamFilters.cpp")
setup_custom_test_program(test_StreamFilters "${SRCROOT}${INPUTOUTPUTDIR}")
target_link_libraries(test_StreamFilters tudat_input_output tudat_basic_astrodynamics ${Boost_LIBRARIES})

add_executable(test_TwoLineElementsTextFileReader "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestTwoLineElementsTextFileReader.cpp")
setup_custom_test_program(test_TwoLineElementsTextFileReader "${SRCROOT}${INPUTOUTPUTDIR}")
target_link_libraries(test_TwoLineElementsTextFileReader tudat_input_output tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_BasicInputOutput "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestBasicInputOutput.cpp")
setup_custom_test_program(test_BasicInputOutput "${SRCROOT}${INPUTOUTPUTDIR}")
target_link_libraries(test_BasicInputOutput tudat_input_output ${Boost_LIBRARIES})

add_executable(test_ParsedDataVectorUtilities "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestParsedDataVectorUtilities.cpp")
setup_custom_test_program(test_ParsedDataVectorUtilities "${SRCROOT}${INPUTOUTPUTDIR}")
target_link_libraries(test_ParsedDataVectorUtilities tudat_input_output ${Boost_LIBRARIES})

add_executable(test_FieldValue "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestFieldValue.cpp")
setup_custom_test_program(test_FieldValue "${SRCROOT}${INPUTOUTPUTDIR}")
target_link_libraries(test_FieldValue tudat_input_output ${Boost_LIBRARIES})

add_executable(test_TextParser "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestTextParser.cpp")
setup_custom_test_program(test_TextParser "${SRCROOT}${INPUTOUTPUTDIR}")
target_link_libraries(test_TextParser tudat_input_output ${Boost_LIBRARIES})

add_executable(test_MissileDatcomReader "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestMissileDatcomReader.cpp")
setup_custom_test_program(test_MissileDatcomReader "${SRCROOT}${INPUTDIR}")
target_link_libraries(test_MissileDatcomReader tudat_input_output ${Boost_LIBRARIES})

add_executable(test_MissileDatcomData "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestMissileDatcomData.cpp")
setup_custom_test_program(test_MissileDatcomData "${SRCROOT}${INPUTDIR}")
target_link_libraries(test_MissileDatcomData tudat_input_output ${Boost_LIBRARIES})

add_executable(test_DictionaryInputSystem "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestDictionaryInputSystem.cpp")
setup_custom_test_program(test_DictionaryInputSystem "${SRCROOT}${INPUTDIR}")
target_link_libraries(test_DictionaryInputSystem tudat_input_output ${Boost_LIBRARIES})

add_executable(test_Extractor "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestExtractor.cpp")
setup_custom_test_program(test_Extractor "${SRCROOT}${INPUTOUTPUTDIR}")
target_link_libraries(test_Extractor tudat_input_output ${Boost_LIBRARIES})

add_executable(test_FixedWidthParser "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestFixedWidthParser.cpp")
setup_custom_test_program(test_FixedWidthParser "${SRCROOT}${INPUTOUTPUTDIR}")
target_link_libraries(test_FixedWidthParser tudat_input_output ${Boost_LIBRARIES})

add_executable(test_SeparatedParser "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestSeparatedParser.cpp")
setup_custom_test_program(test_SeparatedParser "${SRCROOT}${INPUTOUTPUTDIR}")
target_link_libraries(test_SeparatedParser tudat_input_output ${Boost_LIBRARIES})

add_executable(test_LinearFieldTransform "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestLinearFieldTransform.cpp")
setup_custom_test_program(test_LinearFieldTransform "${SRCROOT}${INPUTOUTPUTDIR}")
target_link_libraries(test_LinearFieldTransform tudat_input_output ${Boost_LIBRARIES})

add_executable(test_SolarActivityData "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestSolarActivityData.cpp")
setup_custom_test_program(test_SolarActivityData "${SRCROOT}${INPUTOUTPUTDIR}")
target_link_libraries(test_SolarActivityData tudat_input_output tudat_basic_astrodynamics ${Boost_LIBRARIES})

add_executable(test_MultiArrayReader "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestMultiArrayReader.cpp")
setup_custom_test_program(test_MultiArrayReader "${SRCROOT}${INPUTOUTPUTDIR}")
target_link_libraries(test_MultiArrayReader tudat_input_output tudat_basic_astrodynamics ${Boost_LIBRARIES})

add_executable(test_AerodynamicCoefficientReader "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestAerodynamicCoefficientReader.cpp" )
setup_custom_test_program(test_AerodynamicCoefficientReader "${SRCROOT}${INPUTOUTPUTDIR}")
target_link_libraries(test_AerodynamicCoefficientReader tudat_input_output tudat_basic_astrodynamics tudat_basics ${Boost_LIBRARIES})
