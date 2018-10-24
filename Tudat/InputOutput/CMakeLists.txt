 #    Copyright (c) 2010-2018, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
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
  "${SRCROOT}${INPUTOUTPUTDIR}/streamFilters.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/parseSolarActivityData.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/extractSolarActivityData.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/solarActivityData.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/multiDimensionalArrayReader.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/aerodynamicCoefficientReader.cpp"
  "${SRCROOT}${INPUTOUTPUTDIR}/tabulatedAtmosphereReader.cpp"
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
  "${SRCROOT}${INPUTOUTPUTDIR}/mapTextFileReader.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/matrixTextFileReader.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/streamFilters.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/parseSolarActivityData.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/extractSolarActivityData.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/solarActivityData.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/multiDimensionalArrayReader.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/multiDimensionalArrayWriter.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/aerodynamicCoefficientReader.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/readHistoryFromFile.h"
  "${SRCROOT}${INPUTOUTPUTDIR}/tabulatedAtmosphereReader.h"
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
add_executable(test_MapTextFileReader "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestMapTextFileReader.cpp")
setup_custom_test_program(test_MapTextFileReader "${SRCROOT}${INPUTOUTPUTDIR}")
target_link_libraries(test_MapTextFileReader tudat_input_output tudat_basic_astrodynamics ${Boost_LIBRARIES})

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

add_executable(test_MultiArrayWriter "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestMultiArrayWriter.cpp")
setup_custom_test_program(test_MultiArrayWriter "${SRCROOT}${INPUTOUTPUTDIR}")
target_link_libraries(test_MultiArrayWriter tudat_input_output tudat_basic_astrodynamics ${Boost_LIBRARIES})

add_executable(test_AerodynamicCoefficientReader "${SRCROOT}${INPUTOUTPUTDIR}/UnitTests/unitTestAerodynamicCoefficientReader.cpp" )
setup_custom_test_program(test_AerodynamicCoefficientReader "${SRCROOT}${INPUTOUTPUTDIR}")
target_link_libraries(test_AerodynamicCoefficientReader tudat_input_output tudat_basic_astrodynamics tudat_basics ${Boost_LIBRARIES})
