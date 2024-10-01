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
        "readOdfFile.cpp"
        "readTabulatedMediaCorrections.cpp"
        "readTabulatedWeatherData.cpp"
        "readTrackingTxtFile.cpp"
        "readVariousPdsFiles.cpp"
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
        "readOdfFile.h"
        "readBinaryFile.h"
        "readTabulatedMediaCorrections.h"
        "readTabulatedWeatherData.h"
        "readTrackingTxtFile.h"
        "readVariousPdsFiles.h"
        )

# Add library.
TUDAT_ADD_LIBRARY("input_output"
        "${io_SOURCES}"
        "${io_HEADERS}"
#        PRIVATE_LINKS "${Boost_LIBRARIES}"
#        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )

