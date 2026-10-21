/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/readTrackingTxtFile.h"

#include <cmath>
#include <limits>

namespace tudat
{
namespace input_output
{

static const std::vector< std::string >& getDefaultFdetsDatetimeStringColumnTypes( )
{
    static const std::vector< std::string > columnTypes = {
        "utc_datetime_string", "signal_to_noise_ratio", "normalised_spectral_max", "doppler_measured_frequency_hz", "doppler_noise_hz"
    };
    return columnTypes;
}

static std::vector< std::string > getFdetsColumnTypes( FdetDateFormat dateFormat )
{
    switch( dateFormat )
    {
        case FdetDateFormat::datetime_string:
            return getDefaultFdetsDatetimeStringColumnTypes( );
        case FdetDateFormat::pair_of_numbers:
            throw std::runtime_error( "Fdet files with dates as a pair of values are not currently supported." );
        default:
            throw std::runtime_error( "Invalid Fdet date format." );
    }
}

static std::size_t getFirstFdetsDataLineColumnCount( const std::string& fileName )
{
    std::ifstream dataFile( fileName );
    if( !dataFile.good( ) )
    {
        throw std::runtime_error( "Error when opening Fdets file: file " + fileName + " could not be opened." );
    }

    std::string currentLine;
    while( std::getline( dataFile, currentLine ) )
    {
        boost::algorithm::trim( currentLine );
        if( currentLine.empty( ) || currentLine.at( 0 ) == '#' )
        {
            continue;
        }

        std::vector< std::string > currentSplitRawLine;
        boost::algorithm::split( currentSplitRawLine, currentLine, boost::is_any_of( ", \t" ), boost::algorithm::token_compress_on );
        return currentSplitRawLine.size( );
    }

    return 0;
}

static std::vector< std::string > addOptionalFdetsScanColumnType( const std::string& fileName,
                                                                  const std::vector< std::string >& columnTypes )
{
    if( columnTypes == getDefaultFdetsDatetimeStringColumnTypes( ) && getFirstFdetsDataLineColumnCount( fileName ) == 6 )
    {
        std::vector< std::string > columnTypesWithScanNumber = { "scan_number" };
        columnTypesWithScanNumber.insert( columnTypesWithScanNumber.end( ), columnTypes.begin( ), columnTypes.end( ) );
        return columnTypesWithScanNumber;
    }

    return columnTypes;
}

double getNominalTimeStepFromUtcTimes( const std::vector< double >& observationTimesUtc, const double cadenceTolerance )
{
    if( observationTimesUtc.size( ) < 2 )
    {
        throw std::runtime_error( "Error when getting nominal time step from tracking file contents, size is < 2" );
    }

    double firstObservationTimeStep = std::numeric_limits< double >::infinity( );
    double minimumObservationTimeStep = std::numeric_limits< double >::infinity( );
    for( unsigned int i = 1; i < observationTimesUtc.size( ); i++ )
    {
        double testObservationTimeStep = observationTimesUtc.at( i ) - observationTimesUtc.at( i - 1 );
        if( std::isfinite( testObservationTimeStep ) && testObservationTimeStep > cadenceTolerance )
        {
            if( !std::isfinite( firstObservationTimeStep ) )
            {
                firstObservationTimeStep = testObservationTimeStep;
            }
            if( testObservationTimeStep < minimumObservationTimeStep )
            {
                minimumObservationTimeStep = testObservationTimeStep;
            }
        }
    }

    if( !std::isfinite( minimumObservationTimeStep ) )
    {
        throw std::runtime_error(
                "Error when getting nominal time step from tracking file contents, no positive cadence could be inferred" );
    }

    if( firstObservationTimeStep > minimumObservationTimeStep + cadenceTolerance )
    {
        return minimumObservationTimeStep;
    }
    else
    {
        return firstObservationTimeStep;
    }
}

void TrackingTxtFileContents::parseData( const TrackingTxtFileReadFilterType dataFilterMethod )
{
    std::ifstream dataFile( fileName_ );
    if( !dataFile.good( ) )
    {
        throw std::runtime_error( "Error when opening Tracking txt file: file " + fileName_ + " could not be opened." );
    }
    readRawDataMap( dataFile, dataFilterMethod );
    convertDataMap( dataFilterMethod );
}

void TrackingTxtFileContents::readRawDataMap( std::ifstream& dataFile, const TrackingTxtFileReadFilterType dataFilterMethod )
{
    std::string currentLine;

    while( std::getline( dataFile, currentLine ) )
    {
        if( !currentLine.empty( ) && currentLine.at( 0 ) != commentSymbol_ )
        {
            addLineToRawDataMap( currentLine, dataFilterMethod );
        }
    }
}

bool isIfmsEntryValid( const std::string ifmsFileEntry )
{
    if( ifmsFileEntry.size( ) == 0 )
    {
        return false;
    }
    else
    {
        if( ifmsFileEntry.at( 0 ) != '-' )
        {
            return true;
        }
        else
        {
            bool seenNonNine = false;
            bool seenDot = false;
            for( unsigned int i = 1; i < ifmsFileEntry.size( ); i++ )
            {
                if( ifmsFileEntry.at( i ) == '.' )
                {
                    seenDot = true;
                }
                else if( ifmsFileEntry.at( i ) != '9' )
                {
                    seenNonNine = true;
                }
            }
            return !( seenDot == true && seenNonNine == false );
        }
    }
}

bool TrackingTxtFileContents::validateCurrentLineProcessing( const TrackingTxtFileReadFilterType dataFilterMethod,
                                                             const std::vector< std::string >& rawVector )
{
    bool addLine = true;

    if( dataFilterMethod == ifms_tracking_txt_file_filter )
    {
        for( unsigned int i = 0; i < rawVector.size( ); i++ )
        {
            if( ( i == 6 ) || ( i == 8 ) || ( i == 10 ) )
            {
                addLine = isIfmsEntryValid( rawVector.at( i ) );
                if( !addLine )
                {
                    return addLine;
                }
            }
        }
    }

    return addLine;
}

void TrackingTxtFileContents::addLineToRawDataMap( std::string& rawLine, const TrackingTxtFileReadFilterType dataFilterMethod )
{
    size_t numberOfColumns = getNumColumns( );
    std::vector< std::string > currentSplitRawLine_;

    // Trim the line and split based on the separators
    boost::algorithm::trim( rawLine );
    boost::algorithm::split( currentSplitRawLine_, rawLine, boost::is_any_of( valueSeparators_ ), boost::algorithm::token_compress_on );

    // Check if the expected number of columns is present in this line
    if( currentSplitRawLine_.size( ) != numberOfColumns )
    {
        if( !( currentSplitRawLine_.size( ) > numberOfColumns && ignoreOmittedColumns_ ) )
        {
            unsigned int columnsFound = currentSplitRawLine_.size( );
            for( auto a : currentSplitRawLine_ )
            {
                std::cout << a << "\n";
            }
            throw std::runtime_error( "The current line in file " + fileName_ + " has " + std::to_string( columnsFound ) + " columns but " +
                                      std::to_string( numberOfColumns ) + " columns were expected.\nRaw line:" + rawLine );
        }
    }

    // Populate the dataMap_ with a new row on each of the vectors
    if( validateCurrentLineProcessing( dataFilterMethod, currentSplitRawLine_ ) )
    {
        for( std::size_t i = 0; i < numberOfColumns; ++i )
        {
            std::string currentFieldType = columnFieldTypes_.at( i );
            std::string currentValue = currentSplitRawLine_.at( i );
            rawDataMap_[ currentFieldType ].push_back( currentValue );
        }
    }
}

void TrackingTxtFileContents::convertDataMap( const TrackingTxtFileReadFilterType dataFilterMethod )
{
    // Loop over all the column types to convert them to doubles
    for( std::string columnType : columnFieldTypes_ )
    {
        // If the columnType (requested by the user) does not have a known converter in tudat, it is skipped here and only stored as raw
        // string
        if( !trackingFileFieldConverterMap.count( columnType ) )
        {
            std::cout << "Warning: '" << columnType
                      << "' is not recognised as a column type by Tudat. The data is available in the raw format.\n";
            continue;
        }

        // Get the raw vector, correct converter and target data type represented by the doubles
        std::shared_ptr< TrackingFileFieldConverter > converter = trackingFileFieldConverterMap.at( columnType );
        const std::vector< std::string >& rawVector = rawDataMap_.at( columnType );
        const TrackingDataType& dataType = converter->getTrackingDataType( );

        // Convert and store double vector
        std::vector< double > dataVector;
        for( std::string rawValue : rawVector )
        {
            dataVector.push_back( converter->toDouble( rawValue ) );
        }
        doubleDataMap_[ dataType ] = dataVector;
    }
}

const std::vector< double > TrackingTxtFileContents::getDoubleDataColumn( TrackingDataType dataType, double defaultVal )
{
    if( doubleDataMap_.find( dataType ) != doubleDataMap_.end( ) )
    {
        return doubleDataMap_.at( dataType );
    }
    else if( metaDataMapDouble_.find( dataType ) != metaDataMapDouble_.end( ) )
    {
        return std::vector< double >( getNumRows( ), metaDataMapDouble_[ dataType ] );
    }
    else if( !std::isnan( defaultVal ) )
    {
        return std::vector< double >( getNumRows( ), defaultVal );
    }
    else
    {
        throw std::runtime_error( "Error while working with trackinig data. Required datatype not found" +
                                  std::to_string( static_cast< int >( dataType ) ) );
    }
}

const std::vector< double > TrackingTxtFileContents::getDoubleDataColumn( TrackingDataType dataType )
{
    // Call the version with default value, using a default value that triggers an error
    return getDoubleDataColumn( dataType, std::numeric_limits< double >::quiet_NaN( ) );
}

const std::vector< TrackingDataType > TrackingTxtFileContents::getMetaDataTypes( )
{
    std::vector< TrackingDataType > metaDataTypes;

    for( const auto& pair : metaDataMapDouble_ )
    {
        metaDataTypes.push_back( pair.first );
    }

    for( const auto& pair : metaDataMapStr_ )
    {
        metaDataTypes.push_back( pair.first );
    }

    return metaDataTypes;
}

const std::vector< TrackingDataType > TrackingTxtFileContents::getAllAvailableDataTypes( )
{
    auto columnTypes = getDataColumnTypes( );
    auto metaTypes = getMetaDataTypes( );
    columnTypes.insert( columnTypes.end( ), metaTypes.begin( ), metaTypes.end( ) );
    return columnTypes;
}

void TrackingTxtFileContents::subtractColumnType( const TrackingDataType& columnToSubtractFrom, const TrackingDataType& columnToSubtract )
{
    if( doubleDataMap_.count( columnToSubtractFrom ) == 0 )
    {
        throw std::runtime_error( "Error when subtracing data for tracking text file, column type " +
                                  std::to_string( static_cast< int >( columnToSubtractFrom ) ) + " not found." );
    }

    if( doubleDataMap_.count( columnToSubtract ) == 0 )
    {
        throw std::runtime_error( "Error when subtracing data for tracking text file, column type " +
                                  std::to_string( static_cast< int >( columnToSubtract ) ) + " not found." );
    }

    if( doubleDataMap_.at( columnToSubtractFrom ).size( ) != doubleDataMap_.at( columnToSubtract ).size( ) )
    {
        throw std::runtime_error( "Error when subtracing data for tracking text file, column types " +
                                  std::to_string( static_cast< int >( columnToSubtract ) ) + " and " +
                                  std::to_string( static_cast< int >( columnToSubtractFrom ) ) + " do not have same size." );
    }

    for( unsigned int i = 0; i < doubleDataMap_.at( columnToSubtractFrom ).size( ); i++ )
    {
        doubleDataMap_[ columnToSubtractFrom ][ i ] -= doubleDataMap_[ columnToSubtract ][ i ];
    }
}

std::shared_ptr< TrackingTxtFileContents > readFdetsFile( const std::string& fileName, FdetDateFormat dateFormat )
{
    const std::vector< std::string > columnTypes = addOptionalFdetsScanColumnType( fileName, getFdetsColumnTypes( dateFormat ) );
    auto rawFileContents = createTrackingTxtFileContents( fileName, columnTypes, '#', ", \t" );
    rawFileContents->addMetaData( TrackingDataType::file_name, fileName );
    return rawFileContents;
}

std::shared_ptr< TrackingTxtFileContents > readFdetsFile( const std::string& fileName, const std::vector< std::string >& columnTypes )
{
    static bool hasPrintedDeprecationWarning = false;
    if( !hasPrintedDeprecationWarning )
    {
        std::cerr << "Warning: readFdetsFile(fileName, columnTypes) is deprecated. Use readFdetsFile(fileName, FdetDateFormat) "
                     "instead."
                  << std::endl;
        hasPrintedDeprecationWarning = true;
    }

    const std::vector< std::string > resolvedColumnTypes = addOptionalFdetsScanColumnType( fileName, columnTypes );
    auto rawFileContents = createTrackingTxtFileContents( fileName, resolvedColumnTypes, '#', ", \t" );
    rawFileContents->addMetaData( TrackingDataType::file_name, fileName );
    return rawFileContents;
}

}  // namespace input_output

}  // namespace tudat
