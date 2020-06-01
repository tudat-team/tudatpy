#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# Add source files.
set(io_SOURCES
        "basicInputOutput.cpp"
        "dictionaryComparer.cpp"
        "dictionaryTools.cpp"
        "fieldValue.cpp"
        "fixedWidthParser.cpp"
        "linearFieldTransform.cpp"
        "missileDatcomData.cpp"
        "missileDatcomReader.cpp"
        "parsedDataVectorUtilities.cpp"
        "separatedParser.cpp"
        "textParser.cpp"
        "twoLineElementData.cpp"
        "twoLineElementsTextFileReader.cpp"
        "streamFilters.cpp"
        "parseSolarActivityData.cpp"
        "extractSolarActivityData.cpp"
        "solarActivityData.cpp"
        "multiDimensionalArrayReader.cpp"
        "aerodynamicCoefficientReader.cpp"
        "tabulatedAtmosphereReader.cpp"
        "util.cpp"
        )

# Add header files.
set(io_HEADERS
        "basicInputOutput.h"
        "dictionaryComparer.h"
        "dictionaryEntry.h"
        "dictionaryTools.h"
        "extractor.h"
        "fieldTransform.h"
        "fieldType.h"
        "fieldValue.h"
        "fixedWidthParser.h"
        "linearFieldTransform.h"
        "missileDatcomData.h"
        "missileDatcomReader.h"
        "parsedDataVectorUtilities.h"
        "parser.h"
        "separatedParser.h"
        "textParser.h"
        "twoLineElementData.h"
        "twoLineElementsTextFileReader.h"
        "basicInputOutput.h"
        "mapTextFileReader.h"
        "matrixTextFileReader.h"
        "streamFilters.h"
        "parseSolarActivityData.h"
        "extractSolarActivityData.h"
        "solarActivityData.h"
        "multiDimensionalArrayReader.h"
        "multiDimensionalArrayWriter.h"
        "aerodynamicCoefficientReader.h"
        "readHistoryFromFile.h"
        "tabulatedAtmosphereReader.h"
        "util.h"
        )

# Add library.
TUDAT_ADD_LIBRARY("input_output"
        "${io_SOURCES}"
        "${io_HEADERS}"
        PRIVATE_LINKS "${Boost_LIBRARIES}"
        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )

# Add unit test files.
#set(io_UNITTESTS
#  "${SRCROOT}${INPUTOUTPUTDIR}/tests/unitTestBasicInputOutput.cpp"
#  "${SRCROOT}${INPUTOUTPUTDIR}/tests/unitTestMatrixTextFileReader.cpp"
#  "${SRCROOT}${INPUTOUTPUTDIR}/tests/unitTestStreamFilters.cpp"
#)

TUDAT_ADD_TEST_CASE(MapTextFileReader
        PRIVATE_LINKS
        tudat_input_output
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(MatrixTextFileReader
        PRIVATE_LINKS
        tudat_input_output
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(StreamFilters
        PRIVATE_LINKS
        tudat_input_output
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(TwoLineElementsTextFileReader
        PRIVATE_LINKS
        tudat_input_output
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(BasicInputOutput
        PRIVATE_LINKS
        tudat_input_output
        )

TUDAT_ADD_TEST_CASE(ParsedDataVectorUtilities
        PRIVATE_LINKS
        tudat_input_output
        )

TUDAT_ADD_TEST_CASE(FieldValue
        PRIVATE_LINKS
        tudat_input_output
        )

TUDAT_ADD_TEST_CASE(TextParser
        PRIVATE_LINKS
        tudat_input_output
        )

TUDAT_ADD_TEST_CASE(MissileDatcomReader
        PRIVATE_LINKS
        tudat_input_output
        )

TUDAT_ADD_TEST_CASE(MissileDatcomData
        PRIVATE_LINKS
        tudat_input_output
        )

TUDAT_ADD_TEST_CASE(DictionaryInputSystem
        PRIVATE_LINKS
        tudat_input_output
        )

TUDAT_ADD_TEST_CASE(Extractor
        PRIVATE_LINKS
        tudat_input_output
        )

TUDAT_ADD_TEST_CASE(FixedWidthParser
        PRIVATE_LINKS
        tudat_input_output
        )

TUDAT_ADD_TEST_CASE(SeparatedParser
        PRIVATE_LINKS
        tudat_input_output
        )

TUDAT_ADD_TEST_CASE(LinearFieldTransform
        PRIVATE_LINKS
        tudat_input_output
        )

TUDAT_ADD_TEST_CASE(SolarActivityData
        PRIVATE_LINKS
        tudat_input_output
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(MultiArrayReader
        PRIVATE_LINKS
        tudat_input_output
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(MultiArrayWriter
        PRIVATE_LINKS
        tudat_input_output
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(AerodynamicCoefficientReader
        PRIVATE_LINKS
        tudat_input_output
        tudat_basic_astrodynamics
        tudat_basics
        )
