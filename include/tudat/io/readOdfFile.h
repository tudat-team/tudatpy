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
#include "tudat/io/readBinaryFile.h"
#include "tudat/math/interpolators/lookupScheme.h"

namespace tudat
{
namespace input_output
{

// TODO: test
class OdfClockOffsetBlock
{
public:
    OdfClockOffsetBlock( const std::bitset< 288 > dataBits );

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

    int getPrimaryStationId_( )
    {
        return primaryStationId_;
    }

    int getSecondaryStationId_( )
    {
        return secondaryStationId_;
    }

private:

    unsigned int integerStartTime_; // sec
    unsigned int fractionalStartTime_; // nsec

    int integerClockOffset_;  // sec
    int fractionalClockOffset_; // nsec

    unsigned int primaryStationId_;
    unsigned int secondaryStationId_;

    unsigned int reservedBlock_;

    unsigned int integerEndTime_; // sec
    unsigned int fractionalEndTime_; // nsec
};

class OdfRampBlock
{
public:
    OdfRampBlock( const std::bitset< 288 > dataBits );

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

    int getTransmittingStationId( )
    {
        return transmittingStationId_;
    }

private:

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

};

class OdfDataSpecificBlock
{
public:
    OdfDataSpecificBlock( int dataType_ ): dataType_( dataType_ ){ }

    virtual ~OdfDataSpecificBlock( ){ }

    int dataType_;
};

// TODO: test
// Delta differential one-way Doppler data
class OdfDDodDataBlock: public OdfDataSpecificBlock
{
public:
    OdfDDodDataBlock( const std::bitset< 128 > specificDataBits, const int dDodDataType );

    ~OdfDDodDataBlock( ){ }

    int getSecondReceivingStationId( )
    {
        return secondReceivingStationId_;
    }

    int getQuasarOrSpacecraftId_( )
    {
        return quasarOrSpacecraftId_;
    }

    int getPhasePointIndicator( )
    {
        return phasePointIndicator_;
    }

    double getReferenceFrequency( )
    {
        return std::pow( 2.0, 24 )  / 1.0E3 * referenceFrequencyHighPart_ + referenceFrequencyLowPart_ / 1.0E3;
    }

    int getComposite1( )
    {
        return composite1_;
    }

    double getCompressionTime( )
    {
        return compressionTime_ * 1.0e-2;
    }

    double getSecondReceivingStationUplinkDelay( )
    {
        return secondReceivingStationDownlinkDelay_ * 1.0e-9;
    }

private:

    int secondReceivingStationId_;
    int quasarOrSpacecraftId_;
    int phasePointIndicator_;
    int referenceFrequencyHighPart_; // 2^24 mHz
    int referenceFrequencyLowPart_; // mHz
    int composite1_;
    int compressionTime_; // 1e-2 sec
    int secondReceivingStationDownlinkDelay_; // nsec
};

// TODO: test
// Delta differential one-way ranging data
class OdfDDorDataBlock: public OdfDataSpecificBlock
{
public:
    OdfDDorDataBlock( const std::bitset< 128 > specificDataBits, const int dDorDataType );

    ~OdfDDorDataBlock( ){ }

    int getSecondReceivingStationId( )
    {
        return secondReceivingStationId_;
    }

    int getQuasarOrSpacecraftId_( )
    {
        return quasarOrSpacecraftId_;
    }

    int getModulusIndicator( )
    {
        return modulusIndicator_;
    }

    double getReferenceFrequency( )
    {
        return std::pow( 2.0, 24 )  / 1.0E3 * referenceFrequencyHighPart_ + referenceFrequencyLowPart_ / 1.0E3;
    }

    int getComposite1( )
    {
        return composite1_;
    }

    double getSecondReceivingStationUplinkDelay( )
    {
        return secondReceivingStationDownlinkDelay_ * 1.0e-9;
    }

private:

    int secondReceivingStationId_;
    int quasarOrSpacecraftId_;
    int modulusIndicator_;
    int referenceFrequencyHighPart_; // 2^24 mHz
    int referenceFrequencyLowPart_; // mHz
    int composite1_;
    int modulusLowPart_; // 1e-7 sec
    int secondReceivingStationDownlinkDelay_; // nsec
};

class OdfDopplerDataBlock: public OdfDataSpecificBlock
{
public:
    OdfDopplerDataBlock( const std::bitset< 128 > specificDataBits, const int dopplerDataType );

    ~OdfDopplerDataBlock( ){ }

    int getReceiverChannel( )
    {
        return receiverChannel_;
    }

    int getSpacecraftId( )
    {
        return spacecraftId_;
    }

    int getReceiverExciterFlag( )
    {
        return receiverExciterFlag_;
    }

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

private:

    int receiverChannel_;
    int spacecraftId_;
    int receiverExciterFlag_;
    int referenceFrequencyHighPart_; // 2^24 mHz
    int referenceFrequencyLowPart_; // mHz

    int reservedSegment_;
    int compressionTime_; // 1e-2 sec

    int transmittingStationUplinkDelay_; // nsec

};

class OdfSequentialRangeDataBlock: public OdfDataSpecificBlock
{
public:

    OdfSequentialRangeDataBlock( const std::bitset< 128 > dataBits );

    ~OdfSequentialRangeDataBlock( ){ }

    int getSpacecraftId( )
    {
        return spacecraftId_;
    }

    double getReferenceFrequency( )
    {
        return std::pow( 2.0, 24 ) / 1.0E3 * referenceFrequencyHighPart_ + referenceFrequencyLowPart_ / 1.0E3;
    }

    double getTransmittingStationUplinkDelay( )
    {
        return transmittingStationUplinkDelay_ * 1.0e-9;
    }

    // TODO: create getters for these... need to read about what to get from them though
    int lowestRangingComponent_;
    int reservedBlock_;
    int uplinkCoderInPhaseTimeOffset_; // sec
    int compositeTwo_; // sec

private:

    int spacecraftId_;

    int referenceFrequencyHighPart_; // 2^24 mHz
    int referenceFrequencyLowPart_; // mHz

    int transmittingStationUplinkDelay_; // nsec

};

// TODO: test
class OdfToneRangeDataBlock: public OdfDataSpecificBlock
{
public:

    OdfToneRangeDataBlock( const std::bitset< 128 > specificDataBits );

    ~OdfToneRangeDataBlock( ){ }

    int getSpacecraftId( )
    {
        return spacecraftId_;
    }

    double getReferenceFrequency( )
    {
        return std::pow( 2.0, 24 ) / 1.0E3 * referenceFrequencyHighPart_ + referenceFrequencyLowPart_ / 1.0E3;
    }

    double getTransmittingStationUplinkDelay( )
    {
        return transmittingStationUplinkDelay_ * 1.0e-9;
    }

private:

    int integerObservableTime_; // sec
    int spacecraftId_;
    int reservedBlock1_;

    int referenceFrequencyHighPart_; // 2^24 mHz
    int referenceFrequencyLowPart_; // mHz

    int reservedBlock2_;
    int reservedBlock3_;

    int transmittingStationUplinkDelay_; // nsec
};

// TODO: test
class OdfAngleDataBlock: public OdfDataSpecificBlock
{
public:

    OdfAngleDataBlock( const std::bitset< 128 > specificDataBits, const int angleDataType );

    ~OdfAngleDataBlock( ){ }

    int getSpacecraftId( )
    {
        return spacecraftId_;
    }

private:

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
    OdfCommonDataBlock( const std::bitset< 160 > commonDataBits );

    double getObservableTime( )
    {
        return static_cast< double >( integerTimeTag_ ) + static_cast< double >( fractionalTimeTag_ ) / 1000.0;
    }

    double getObservableValue( )
    {
        return static_cast< double >( integerObservable_ ) + static_cast< double >( fractionalObservable_ ) / 1.0E9;
    }

    double getReceivingStationDownlinkDelay( )
    {
        return receivingStationDownlinkDelay_ * 1.0e-9;
    }

    int formatId_;
    int receivingStationId_;
    int transmittingStationId_;
    int transmittingStationNetworkId_;
    int dataType_;
    int downlinkBandId_;
    int uplinkBandId_;
    int referenceBandId_;
    int validity_;

private:

    uint32_t integerTimeTag_; // sec
    int fractionalTimeTag_; // msec
    int receivingStationDownlinkDelay_; // nsec

    int integerObservable_; // unit
    int fractionalObservable_; // 1e-9 * unit
};

class OdfDataBlock
{
public:
    OdfDataBlock( const std::bitset< 288 > dataBits );

    std::shared_ptr< OdfDataSpecificBlock > getObservableSpecificDataBlock( )
    {
        return observableSpecificDataBlock_;
    }

    std::shared_ptr< OdfCommonDataBlock > getCommonDataBlock( )
    {
        return commonDataBlock_;
    }

private:

    std::shared_ptr< OdfDataSpecificBlock > observableSpecificDataBlock_;
    std::shared_ptr< OdfCommonDataBlock > commonDataBlock_;

};

class OdfRawFileContents
{
public:
    /*!
     *
     * @param odfFile File name/location of ODF file that is to be read
     */
    OdfRawFileContents( const std::string& odfFile );

    std::string systemId_;
    std::string programId_;
    uint32_t spacecraftId_;

    uint32_t fileCreationDate_; // year, month, day (YYYMMDD): year from 1900 or 1950
    uint32_t fileCreationTime_; // hour, minute, second (HHMMSS)

    uint32_t fileReferenceDate_; // year, month, day (YYYYMMDD)
    uint32_t fileReferenceTime_; // hour, minute, second (HHMMSS)

    std::string fileName_;

    std::string identifierGroupStringA_;
    std::string identifierGroupStringB_;
    std::string identifierGroupStringC_;

    bool eofHeaderFound_;

    //! Function to retrieve the orbit data blocsk
    std::vector< std::shared_ptr< OdfDataBlock > > getDataBlocks( )
    {
        return dataBlocks_;
    }

    //! Function to retrieve the ramp blocks
    std::map< int, std::vector< std::shared_ptr< OdfRampBlock > > > getRampBlocks( )
    {
        return rampBlocks_;
    }

    //! Function to retrieve the clock offset blocks
    std::map< std::pair< int, int >, std::shared_ptr< OdfClockOffsetBlock > > getClockOffsetBlocks( )
    {
        return clockOffsetBlocks_;
    }

private:

    //! Vector of data blocks
    std::vector< std::shared_ptr< OdfDataBlock > > dataBlocks_;

    //! Vector of ramp blocks indexed by transmitting station ID
    std::map< int, std::vector< std::shared_ptr< OdfRampBlock > > > rampBlocks_;

    //! Clock offset blocks indexed by pair of (primary station ID, secondary station ID)
    std::map< std::pair< int, int >, std::shared_ptr< OdfClockOffsetBlock > > clockOffsetBlocks_;

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

};

} // namespace input_output

} // namespace tudat

#endif // TUDAT_READ_ODF_FILE_H
