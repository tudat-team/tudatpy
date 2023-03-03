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

    unsigned int integerRampStartTime_; // sec
    unsigned int fractionalRampStartTime_; // nsec

    int integerRampRate_; // Hz/s
    int fractionalRampRate_; // 1e-9 Hz/s

    int integerRampStartFrequency_;  // GHz

    int transmittingStationId_;

    unsigned int integerRampStartFrequencyModulo_; // Hz
    unsigned int fractionalRampStartFrequency_; // Hz

    unsigned int integerRampEndTime_; // sec
    unsigned int fractionalRampEndTime_; // nsec

    double getRampStartFrequency( )
    {
        return static_cast< double >( integerRampStartFrequency_ ) * 1.0E9 +
               static_cast< double >( integerRampStartFrequencyModulo_ ) +
               static_cast< double >( fractionalRampStartFrequency_ ) * 1.0E-9;
    }

    double getRampRate( )
    {
        return static_cast< double >( integerRampRate_ ) +
               static_cast< double >( fractionalRampRate_ ) * 1.0E-9;
    }

    double getRampStartTime( )
    {
        return static_cast< double >( integerRampStartTime_ ) + static_cast< double >( fractionalRampStartTime_ ) * 1.0E-9;
    }

    double getRampEndTime( )
    {
        return static_cast< double >( integerRampEndTime_ ) + static_cast< double >( fractionalRampEndTime_ ) * 1.0E-9;
    }

    void printContents( )
    {
        std::cout << "Start time " << integerRampStartTime_ << " " << fractionalRampStartTime_ << std::endl;
        std::cout << "End time " << integerRampEndTime_ << " " << fractionalRampEndTime_ << std::endl;
        std::cout << "Ramp rate " << integerRampRate_ << " " << fractionalRampRate_ << std::endl;
        std::cout << "Start frequency " << integerRampStartFrequency_ << " " << integerRampStartFrequencyModulo_ << " " <<
                  fractionalRampStartFrequency_ << std::endl;
        std::cout << "Station " << transmittingStationId_ << std::endl << std::endl;
    }
};

class OdfDataSpecificBlock
{
public:
    OdfDataSpecificBlock( int dataType_ ): dataType_( dataType_ ){ }

    virtual ~OdfDataSpecificBlock( ){ }

    int dataType_;
};

// Delta differential one-way Doppler data
class OdfDDodDataBlock: public OdfDataSpecificBlock
{
public:
    OdfDDodDataBlock( const int dDodDataType ):
        OdfDataSpecificBlock( dDodDataType ){ }

    ~OdfDDodDataBlock( ){ }

    int secondReceivingStationId_;
    int quasarOrSpacecraftId_;
    int phasePointIndicator_;
    int referenceFrequencyHighPart_; // 2^24 mHz
    int referenceFrequencyLowPart_; // mHz
    int composite1_;
    int compressionTime_; // 1e-2 sec
    int secondReceivingStationDownlinkDelay_; // nsec

    double getReferenceFrequency( )
    {
        return std::pow( 2.0, 24 )  / 1.0E3 * referenceFrequencyHighPart_ + referenceFrequencyLowPart_ / 1.0E3;
    }
};

// Delta differential one-way ranging data
class OdfDDorDataBlock: public OdfDataSpecificBlock
{
public:
    OdfDDorDataBlock( const int dDorDataType ):
        OdfDataSpecificBlock( dDorDataType ){ }

    ~OdfDDorDataBlock( ){ }

    int secondReceivingStationId_;
    int quasarOrSpacecraftId_;
    int modulusIndicator_;
    int referenceFrequencyHighPart_; // 2^24 mHz
    int referenceFrequencyLowPart_; // mHz
    int composite1_;
    int modulusLowPart_; // 1e-7 sec
    int secondReceivingStationDownlinkDelay_; // nsec

    double getReferenceFrequency( )
    {
        return std::pow( 2.0, 24 )  / 1.0E3 * referenceFrequencyHighPart_ + referenceFrequencyLowPart_ / 1.0E3;
    }
};

class OdfDopplerDataBlock: public OdfDataSpecificBlock
{
public:
    OdfDopplerDataBlock( const int DopplerDataType ):
        OdfDataSpecificBlock( DopplerDataType ){ }

    ~OdfDopplerDataBlock( ){ }

    int receiverChannel_;
    int spacecraftId_;
    int receiverExciterFlag_;
    int referenceFrequencyHighPart_; // 2^24 mHz
    int referenceFrequencyLowPart_; // mHz

    int reservedSegment_;
    int compressionTime_; // 1e-2 sec

    int transmittingStationUplinkDelay_; // nsec

    double getReferenceFrequency( )
    {
        return std::pow( 2.0, 24 ) / 1.0E3 * referenceFrequencyHighPart_ + referenceFrequencyLowPart_ / 1.0E3;
    }

    double getCompressionTime( )
    {
        return compressionTime_ * 1.0e-2;
    }

    double getTransmittingStationUplinkDelay( )
    {
        return transmittingStationUplinkDelay_ * 1.0e-9;
    }

    void printContents( )
    {
        std::cout << "Receiver: " << receiverChannel_ << " " << receiverExciterFlag_ << std::endl;
        std::cout << "Spacecraft: " << spacecraftId_ << std::endl;
        std::cout << "Reference frequency: " << referenceFrequencyHighPart_ << " " << referenceFrequencyLowPart_ << std::endl;

        std::cout << "Reserved: " << reservedSegment_ << std::endl;
        std::cout << "Compression time: " << compressionTime_ << std::endl;
        std::cout << "Transmission delay: " << transmittingStationUplinkDelay_ << std::endl << std::endl;;

    }
};

class OdfSequentialRangeDataBlock: public OdfDataSpecificBlock
{
public:

    OdfSequentialRangeDataBlock( ):
        OdfDataSpecificBlock( 37 ){ }

    ~OdfSequentialRangeDataBlock( ){ }

    int lowestRangingComponent_;
    int spacecraftId_;
    int reservedBlock_;

    int referenceFrequencyHighPart_; // 2^24 mHz
    int referenceFrequencyLowPart_; // mHz

    int coderInPhaseTimeOffset_; // sec
    int compositeTwo_; // sec
    int transmittingStationUplinkDelay_; // nsec

    double getReferenceFrequency( )
    {
        return std::pow( 2.0, 24 ) / 1.0E3 * referenceFrequencyHighPart_ + referenceFrequencyLowPart_ / 1.0E3;
    }

    double getTransmittingStationUplinkDelay( )
    {
        return transmittingStationUplinkDelay_ * 1.0e-9;
    }
};

class OdfToneRangeDataBlock: public OdfDataSpecificBlock
{
public:

    OdfToneRangeDataBlock( ):
        OdfDataSpecificBlock( 41 ){ }

    ~OdfToneRangeDataBlock( ){ }

    int integerObservableTime_; // sec
    int spacecraftId_;
    int reservedBlock1_;

    int referenceFrequencyHighPart_; // 2^24 mHz
    int referenceFrequencyLowPart_; // mHz

    int reservedBlock2_;
    int reservedBlock3_;

    int transmittingStationUplinkDelay_; // nsec

    double getReferenceFrequency( )
    {
        return std::pow( 2.0, 24 ) / 1.0E3 * referenceFrequencyHighPart_ + referenceFrequencyLowPart_ / 1.0E3;
    }
};

class OdfAngleDataBlock: public OdfDataSpecificBlock
{
public:

    OdfAngleDataBlock( const int angleDataType ):
        OdfDataSpecificBlock( angleDataType ){ }

    ~OdfAngleDataBlock( ){ }

    int reservedBlock1_;
    int spacecraftId_;
    int reservedBlock2_;
    int reservedBlock3_;
    int reservedBlock4_;
    int reservedBlock5_;
    int reservedBlock6_;
    int reservedBlock7_;
};

class OdfCommonDataBlock
{
public:

    double getObservableTime( )
    {
        return static_cast< double >( integerTimeTag_ ) + static_cast< double >( fractionalTimeTag_ ) / 1000.0;
    }

    double getObservableValue( )
    {
        return static_cast< double >( integerObservable_ ) +
               static_cast< double >( fractionalObservable_ ) / 1.0E9;
    }

    double getReceivingStationDownlinkDelay( )
    {
        return receivingStationDownlinkDelay_ * 1.0e-9;
    }

    uint32_t integerTimeTag_; // sec
    int fractionalTimeTag_; // msec
    int receivingStationDownlinkDelay_; // nsec

    int integerObservable_; // unit
    int fractionalObservable_; // 1e-9 * unit

    int formatId_;
    int receivingStationId_;
    int transmittingStationId_;
    int transmittingStationNetworkId_;

    int downlinkBandId_;
    int uplinkBandId_;
    int referenceBandId_;
    int validity_;

    void printContents( )
    {
        std::cout << "Time: " << integerTimeTag_ << " " << fractionalTimeTag_ << std::endl;
        std::cout << "Downlink delay: " << receivingStationDownlinkDelay_ << std::endl;
        std::cout << "Observable: " << integerObservable_ << " " << fractionalObservable_ << std::endl;

        std::cout << "Format id: " << formatId_ << std::endl;
        std::cout << "Station data: " << receivingStationId_ << " " << transmittingStationId_ << " " << transmittingStationNetworkId_ << " " << std::endl;
        std::cout << "Bands: " << downlinkBandId_ << " " << uplinkBandId_ << " " << referenceBandId_ << " " << std::endl;
        std::cout << "Validity: " << validity_ << std::endl << std::endl;

    }
};

class OdfDataBlock
{
public:
    std::shared_ptr< OdfDataSpecificBlock > observableSpecificDataBlock_;
    std::shared_ptr< OdfCommonDataBlock > commonDataBlock_;
};

class OdfRawFileContents
{
public:
    std::string systemId_;
    std::string programId_;
    uint32_t spacecraftId_;

    uint32_t fileCreationDate_; // year, month, day (YYMMDD)
    uint32_t fileCreationTime_; // hour, minute, second (HHMMSS)

    uint32_t fileReferenceDate_; // year, month, day (YYMMDD)
    uint32_t fileReferenceTime_; // hour, minute, second (HHMMSS)

    std::string fileName_;

    std::string identifierGroupStringA_;
    std::string identifierGroupStringB_;
    std::string identifierGroupStringC_;

    bool eofHeaderFound_;

    std::vector< std::shared_ptr< OdfDataBlock > > dataBlocks_;

    // Indexed by transmitting station ID
    std::map< int, std::vector< std::shared_ptr< OdfRampBlock > > > rampBlocks_;

    // Indexed by pair of (primary station ID, secondary station ID)
    std::map< std::pair< int, int >, std::shared_ptr< OdfClockOffsetBlock > > clockOffsetBlocks_;
};

std::shared_ptr< OdfClockOffsetBlock > parseClockOffsetData( std::bitset< 288 > dataBits );

//! Function to parse the contents of an ODF orbit data block, specific for sequential range data.
std::shared_ptr< OdfSequentialRangeDataBlock > parseSequentialRangeData( std::bitset< 128 > dataBits );

//! Function to parse the contents of an ODF orbit data block, specific for Doppler data.
std::shared_ptr< OdfDopplerDataBlock > parseDopplerOrbitData( std::bitset< 128 > dataBits, const int dopplerType );

std::shared_ptr< OdfDDodDataBlock > parseDDodOrbitData( std::bitset< 128 > dataBits, const int dDodType );

std::shared_ptr< OdfDDorDataBlock > parseDDorOrbitData( std::bitset< 128 > dataBits, const int dDorType );

std::shared_ptr< OdfToneRangeDataBlock > parseToneRangeOrbitData( std::bitset< 128 > dataBits );

std::shared_ptr< OdfAngleDataBlock > parseAngleOrbitData( std::bitset< 128 > dataBits, const int angleType );


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
