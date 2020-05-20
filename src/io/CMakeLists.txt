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
        "${Boost_LIBRARIES}"
        )

# Add unit test files.
#set(io_UNITTESTS
#  "${SRCROOT}${INPUTOUTPUTDIR}/tests/unitTestBasicInputOutput.cpp"
#  "${SRCROOT}${INPUTOUTPUTDIR}/tests/unitTestMatrixTextFileReader.cpp"
#  "${SRCROOT}${INPUTOUTPUTDIR}/tests/unitTestStreamFilters.cpp"
#)

TUDAT_ADD_TESTCASE(MapTextFileReader
        tudat_input_output
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TESTCASE(MatrixTextFileReader
        tudat_input_output
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TESTCASE(StreamFilters
        tudat_input_output
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TESTCASE(TwoLineElementsTextFileReader
        tudat_input_output
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TESTCASE(BasicInputOutput
        tudat_input_output
        )

TUDAT_ADD_TESTCASE(ParsedDataVectorUtilities
        tudat_input_output
        )

TUDAT_ADD_TESTCASE(FieldValue
        tudat_input_output
        )

TUDAT_ADD_TESTCASE(TextParser
        tudat_input_output
        )

TUDAT_ADD_TESTCASE(MissileDatcomReader
        tudat_input_output
        )

TUDAT_ADD_TESTCASE(MissileDatcomData
        tudat_input_output
        )

TUDAT_ADD_TESTCASE(DictionaryInputSystem
        tudat_input_output
        )

TUDAT_ADD_TESTCASE(Extractor
        tudat_input_output
        )

TUDAT_ADD_TESTCASE(FixedWidthParser
        tudat_input_output
        )

TUDAT_ADD_TESTCASE(SeparatedParser
        tudat_input_output
        )

TUDAT_ADD_TESTCASE(LinearFieldTransform
        tudat_input_output
        )

TUDAT_ADD_TESTCASE(SolarActivityData
        tudat_input_output
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TESTCASE(MultiArrayReader
        tudat_input_output
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TESTCASE(MultiArrayWriter
        tudat_input_output
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TESTCASE(AerodynamicCoefficientReader
        tudat_input_output
        tudat_basic_astrodynamics
        tudat_basics
        )
