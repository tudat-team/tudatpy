/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_READ_ODF_FILE_H
#define TUDAT_READ_ODF_FILE_H

#include <iostream>
#include <fstream>
#include <bitset>
#include <cmath>
#include <vector>
#include <map>

#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/math/interpolators/lookupScheme.h"

namespace tudat
{
namespace input_output
{

class OdfClockOffsetBlock
{
public:

    unsigned int integerStartTime_; // sec
    unsigned int fractionalStartTime_; // nsec

    int integerClockOffset_;  // sec
    int fractionalClockOffset_; // nsec

    unsigned int primaryStationId_;
    unsigned int secondaryStationId_;

    unsigned int reservedBlock_;

    unsigned int integerEndTime_; // sec
    unsigned int fractionalEndTime_; // nsec

    double getStartTime( )
    {
        return static_cast< double >( integerStartTime_ ) + static_cast< double >( fractionalStartTime_ ) * 1.0E-9;
    }

    double getEndTime( )
    {
        return static_cast< double >( integerEndTime_ ) + static_cast< double >( fractionalEndTime_ ) * 1.0E-9;
    }

    double getClockOffset( )
    {
        return static_cast< double >( integerClockOffset_ ) + static_cast< double >( fractionalClockOffset_ ) * 1.0E-9;
    }
};

class OdfRampBlock
{
public:

    unsigned int integerRampStartTime; // sec
    unsigned int fractionalRampStartTime; // nsec

    int integerRampRate; // Hz/s
    int fractionalRampRate; // 1e-9 Hz/s

    int integerRampStartFrequency;  // GHz

    int transmittingStationId;

    unsigned int integerRampStartFrequencyModulo; // Hz
    unsigned int fractionalRampStartFrequency; // Hz

    unsigned int integerRampEndTime; // sec
    unsigned int fractionalRampEndTime; // nsec

    double getRampStartFrequency( )
    {
        return static_cast< double >( integerRampStartFrequency ) * 1.0E9 +
                static_cast< double >( integerRampStartFrequencyModulo ) +
                static_cast< double >( fractionalRampStartFrequency ) * 1.0E-9;
    }

    double getRampRate( )
    {
        return static_cast< double >( integerRampRate ) +
                static_cast< double >( fractionalRampRate ) * 1.0E-9;
    }

    double getRampStartTime( )
    {
        return static_cast< double >( integerRampStartTime ) + static_cast< double >( fractionalRampStartTime ) * 1.0E-9;
    }

    double getRampEndTime( )
    {
        return static_cast< double >( integerRampEndTime ) + static_cast< double >( fractionalRampEndTime ) * 1.0E-9;
    }

    void printContents( )
    {
        std::cout<<"Start time "<<integerRampStartTime<<" "<<fractionalRampStartTime<<std::endl;
        std::cout<<"End time "<<integerRampEndTime<<" "<<fractionalRampEndTime<<std::endl;
        std::cout<<"Ramp rate "<<integerRampRate<<" "<<fractionalRampRate<<std::endl;
        std::cout<<"Start frequency "<<integerRampStartFrequency<<" "<<integerRampStartFrequencyModulo<<" "<<
                   fractionalRampStartFrequency<<std::endl;
        std::cout<<"Station "<<transmittingStationId<<std::endl<<std::endl;
    }
};

class OdfDataSpecificBlock
{
public:
    OdfDataSpecificBlock( int dataType_ ): dataType_( dataType_ ){ }

    virtual ~OdfDataSpecificBlock( ){ }

    int dataType_;
};

class OdfSequentialRangeDataBlock: public OdfDataSpecificBlock
{
public:

    OdfSequentialRangeDataBlock( ):
        OdfDataSpecificBlock( 37 ){ }

    ~OdfSequentialRangeDataBlock( ){ }

    int lowestRangingComponent;
    int spacecraftId;
    int reservedBlock;

    int referenceFrequencyHighPart; // 2^24 mHz
    int referenceFrequencyLowPart; // mHz

    int coderInPhaseTimeOffset; // sec
    int compositeTwo; // sec
    int transmittingStationUplinkDelay; // nsec

    double getReferenceFrequency( )
    {
        return std::pow( 2.0, 24 )  / 1.0E3 * referenceFrequencyHighPart + referenceFrequencyLowPart / 1.0E3;
    }

};

class OdfDopplerDataBlock: public OdfDataSpecificBlock
{
public:
    OdfDopplerDataBlock( const int DopplerDataType ):
        OdfDataSpecificBlock( DopplerDataType ){ }

    ~OdfDopplerDataBlock( ){ }

    int receiverChannel;
    int spacecraftId;
    int receiverExciterFlag;
    int referenceFrequencyHighPart; // 2^24 mHz
    int referenceFrequencyLowPart; // mHz

    int reservedSegment;
    int compressionTime; // 1e-2 sec

    int transmittingStationUplinkDelay; // nsec

    double getReferenceFrequency( )
    {
        return std::pow( 2.0, 24 )  / 1.0E3 * referenceFrequencyHighPart + referenceFrequencyLowPart / 1.0E3;
    }

    void printContents( )
    {
        std::cout<<"Receiver: "<<receiverChannel<<" "<<receiverExciterFlag<<std::endl;
        std::cout<<"Spacecraft: "<<spacecraftId<<std::endl;
        std::cout<<"Reference frequency: "<<referenceFrequencyHighPart<<" "<<referenceFrequencyLowPart<<std::endl;

        std::cout<<"Reserved: "<<reservedSegment<<std::endl;
        std::cout<<"Compression time: "<<compressionTime<<std::endl;
        std::cout << "Transmission delay: " << transmittingStationUplinkDelay << std::endl << std::endl;;

    }

};

class OdfCommonDataBlock
{
public:

    double getObservableValue( )
    {
        return static_cast< double >( integerObservable ) +
                static_cast< double >( fractionalObservable ) / 1.0E9;
    }

    double getObservableTime( )
    {
        return static_cast< double >( integerTimeTag ) + static_cast< double >( fractionalTimeTag ) / 1000.0;
    }

    uint32_t integerTimeTag; // sec
    int fractionalTimeTag; // msec
    int receivingStationDownlinkDelay; // nsec

    int integerObservable; // unit
    int fractionalObservable; // 1e-9 * unit

    int formatId;
    int receivingStation;
    int transmittingStation;
    int transmittingStationNetworkId;

    int downlinkBand;
    int uplinkBand;
    int referenceBand;
    int validity;

    void printContents( )
    {
        std::cout<<"Time: "<<integerTimeTag<<" "<<fractionalTimeTag<<std::endl;
        std::cout<<"Downlink delay: "<<receivingStationDownlinkDelay<<std::endl;
        std::cout<<"Observable: "<<integerObservable<<" "<<fractionalObservable<<std::endl;

        std::cout<<"Format id: "<<formatId<<std::endl;
        std::cout<<"Station data: "<<receivingStation<<" "<<transmittingStation<<" "<<transmittingStationNetworkId<<" "<<std::endl;
        std::cout<<"Bands: "<<downlinkBand<<" "<<uplinkBand<<" "<<referenceBand<<" "<<std::endl;
        std::cout<<"Validity: "<<validity<<std::endl<<std::endl;

    }
};

class OdfDataBlock
{
public:
    std::shared_ptr< OdfDataSpecificBlock > observableSpecificDataBlock;
    std::shared_ptr< OdfCommonDataBlock > commonDataBlock;
};

class OdfRawFileContents
{
public:
    std::string systemId;
    std::string programId;
    uint32_t spacecraftId;

    uint32_t fileCreationDate; // year, month, day (YYMMDD)
    uint32_t fileCreationTime; // hour, minute, second (HHMMSS)

    uint32_t fileReferenceDate; // year, month, day (YYMMDD)
    uint32_t fileReferenceTime; // hour, minute, second (HHMMSS)

    std::string fileName;

    std::string identifierGroupStringA;
    std::string identifierGroupStringB;
    std::string identifierGroupStringC;

    bool eofHeaderFound;

    std::vector< std::shared_ptr< OdfDataBlock > > dataBlocks;

    // Indexed by transmitting station ID
    std::map< int, std::vector< std::shared_ptr< OdfRampBlock > > > rampBlocks;

    // Indexed by pair of (primary station ID, secondary station ID)
    std::map< std::pair< int, int >, std::shared_ptr< OdfClockOffsetBlock > > clockOffsetBlocks;
};

std::shared_ptr< OdfClockOffsetBlock > parseClockOffsetData( std::bitset< 288 > dataBits );

//! Function to parse the contents of an ODF orbit data block, specific for sequenctial range data.
std::shared_ptr< OdfSequentialRangeDataBlock > parseSequentialRangeData( std::bitset< 128 > dataBits );

//! Function to parse the contents of an ODF orbit data block, specific for Doppler data.
std::shared_ptr< OdfDopplerDataBlock > parseDopplerOrbitData( std::bitset< 128 > dataBits, const int dopplerType );

//! Function to parse the contents of an ODF orbit data block
std::shared_ptr< OdfDataBlock > parseOrbitData( std::bitset< 288 > dataBits );

//! Function to parse the contents of an ODF ramp data block
std::shared_ptr< OdfRampBlock > parseRampData( std::bitset< 288 > dataBits );

//! Function to parse the contents of an ODF file label block
void parseFileLabelData(
        std::bitset< 288 > dataBits, std::string& systemId, std::string& programId, uint32_t& spacecraftId,
        uint32_t& fileCreationDate, uint32_t& fileCreationTime, uint32_t& fileReferenceDate,
        uint32_t& fileReferenceTime );

void parseIdentifierData(
        std::bitset< 288 > dataBits, std::string& identifierGroupStringA, std::string&identifierGroupStringB,
        std::string& identifierGroupStringC );

//! Function to parse the contents of an ODF file header block
void parseHeader(
        std::bitset< 288 > dataBits, int32_t& primaryKey, uint32_t& secondaryKey, uint32_t& logicalRecordLength,
        uint32_t& groupStartPacketNumber);

//! Function to check if the current ODF data block is a header.
bool currentBlockIsHeader(
        std::bitset< 288 > dataBits, int& primaryKey, unsigned int& secondaryKey,
        unsigned int& logicalRecordLength, unsigned int& groupStartPacketNumber );

//! Function to read a single 36 byte block from ODF file
void readOdfFileBlock( std::istream& file, std::bitset< 36 * 8 >& dataBits );

//! Function to read the contents of an ODF file into an OdfRawFileContents object
/*!
 * Function to read the contents of an ODF file into an OdfRawFileContents object. The OdfRawFileContents object contains the
 * unprocessed contents of teh ODF file, on a line-by-line basis.
 * \param odfFile File name/location of ODF file that is to be read
 * \return OdfRawFileContents object with contents of ODF file.
 */
std::shared_ptr< OdfRawFileContents > readOdfFile(
        const std::string& odfFile );

} // namespace input_output

} // namespace tudat

#endif // TUDAT_READ_ODF_FILE_H
