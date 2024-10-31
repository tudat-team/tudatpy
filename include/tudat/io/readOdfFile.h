/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References: 820-013, TRK-2-18 Tracking System Interfaces Orbit Data File Interface, Revision
 * E, 2008, JPL/DSN
 *
 */

#ifndef TUDAT_READ_ODF_FILE_H
#define TUDAT_READ_ODF_FILE_H

#include <bitset>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/basics/timeType.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/io/readBinaryFile.h"
#include "tudat/math/interpolators/lookupScheme.h"

namespace tudat
{
namespace input_output
{

// TODO: test
//! ODF file clock offset block
class OdfClockOffsetBlock
{
   public:
    /*!
     * Constructor. Parses an ODF clock offset block, according to table 3-6 of TRK-2-18 (2018).
     *
     * @param dataBits Data block of ODF file.
     */
    OdfClockOffsetBlock( const std::bitset< 288 > dataBits );

    // Returns the start time of the clock offset validity in UTC seconds since the reference time
    // specified in the header.
    double getStartTime( )
    {
        return static_cast< double >( integerStartTime_ ) +
                static_cast< double >( fractionalStartTime_ ) * 1.0E-9;
    }

    // Returns the end time of the clock offset validity in UTC seconds since the reference time
    // specified in the header.
    double getEndTime( )
    {
        return static_cast< double >( integerEndTime_ ) +
                static_cast< double >( fractionalEndTime_ ) * 1.0E-9;
    }

    // Returns the clock offset in seconds.
    double getClockOffset( )
    {
        return static_cast< double >( integerClockOffset_ ) +
                static_cast< double >( fractionalClockOffset_ ) * 1.0E-9;
    }

    // Returns the ID of the primary DSN receiving station.
    int getPrimaryStationId_( )
    {
        return primaryStationId_;
    }

    // Returns the ID of the secondary DSN receiving station.
    int getSecondaryStationId_( )
    {
        return secondaryStationId_;
    }

   private:
    unsigned int integerStartTime_;     // sec
    unsigned int fractionalStartTime_;  // nsec

    int integerClockOffset_;     // sec
    int fractionalClockOffset_;  // nsec

    unsigned int primaryStationId_;
    unsigned int secondaryStationId_;

    unsigned int reservedBlock_;

    unsigned int integerEndTime_;     // sec
    unsigned int fractionalEndTime_;  // nsec
};

//! ODF file ramp block
class OdfRampBlock
{
   public:
    /*!
     * Constructor. Parses an ODF ramp block, according to table 3-5 of TRK-2-18 (2018).
     *
     * @param dataBits Data block of ODF file.
     */
    OdfRampBlock( const std::bitset< 288 > dataBits );

    // Returns the ramp start frequency in Hz.
    double getRampStartFrequency( )
    {
        return static_cast< double >( integerRampStartFrequency_ ) * 1.0E9 +
                static_cast< double >( integerRampStartFrequencyModulo_ ) +
                static_cast< double >( fractionalRampStartFrequency_ ) * 1.0E-9;
    }

    // Returns the constant ramp rate in Hz/s.
    double getRampRate( )
    {
        return static_cast< double >( integerRampRate_ ) +
                static_cast< double >( fractionalRampRate_ ) * 1.0E-9;
    }

    // Returns the ramp start time in UTC seconds since the reference time specified in the header.
    Time getRampStartTime( )
    {
        return Time( static_cast< double >( integerRampStartTime_ ) ) +
                Time( static_cast< double >( fractionalRampStartTime_ ) * 1.0E-9 );
    }

    // Returns the ramp end time in UTC seconds since the reference time specified in the header.
    Time getRampEndTime( )
    {
        return Time( static_cast< double >( integerRampEndTime_ ) ) +
                Time( static_cast< double >( fractionalRampEndTime_ ) * 1.0E-9 );
    }

    // ID of the DSN station to which the ramp block applies, specified according to TRK-2-18
    // (2018).
    int getTransmittingStationId( )
    {
        return transmittingStationId_;
    }

    /*!
     * Prints the data block to the specified output file.
     *
     * @param outFile File name on which to write data.
     */
    void printDataBlock( std::ofstream &outFile );

   private:
    unsigned int integerRampStartTime_;     // sec
    unsigned int fractionalRampStartTime_;  // nsec

    int integerRampRate_;     // Hz/s
    int fractionalRampRate_;  // 1e-9 Hz/s

    int integerRampStartFrequency_;  // GHz

    // ID of the DSN station to which the ramp block applies
    int transmittingStationId_;

    unsigned int integerRampStartFrequencyModulo_;  // Hz
    unsigned int fractionalRampStartFrequency_;     // Hz

    unsigned int integerRampEndTime_;     // sec
    unsigned int fractionalRampEndTime_;  // nsec
};

// Base class defining the observable specific portion of an ODF data block.
class OdfDataSpecificBlock
{
   public:
    /*!
     * Constructor.
     * @param dataType_ Data type, specified according to section 3.2.4 of TRK-2-18 (2018).
     */
    OdfDataSpecificBlock( int dataType_ ) : dataType_( dataType_ ) { }

    // Destructor
    virtual ~OdfDataSpecificBlock( ) { }

    /*!
     * Virtual function printing the observable specific data block to the specified output file.
     * Should be implemented in the derived classes.
     *
     * @param outFile File name on which to write data.
     */
    virtual void printDataBlock( std::ofstream &outFile )
    {
        outFile << "Printing not implemented for current type";
    }

    // Observable type, according to table 3-4a, item 10, of TRK-2-18 (2018).
    int dataType_;
};

// TODO: test
// Derived class defining a Delta differential one-way Doppler data block
class OdfDDodDataBlock : public OdfDataSpecificBlock
{
   public:
    /*!
     * Constructor. Parses a the observable specific portion of a Delta differential one-way Doppler
     * block, according to table 3-4b of TRK-2-18 (2018).
     *
     * @param specificDataBits Observable specific portion of the data block.
     * @param dDodDataType Data type, specified according to section 3.2.4 of TRK-2-18 (2018). Can take values 1, 2, 3, 4.
     */
    OdfDDodDataBlock( const std::bitset< 128 > specificDataBits, const int dDodDataType );

    // Destructor
    ~OdfDDodDataBlock( ) { }

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

    // Returns the reference frequency in Hz
    double getReferenceFrequency( )
    {
        return std::pow( 2.0, 24 ) / 1.0E3 * referenceFrequencyHighPart_ +
                referenceFrequencyLowPart_ / 1.0E3;
    }

    int getComposite1( )
    {
        return composite1_;
    }

    // Returns the compression time in seconds.
    double getCompressionTime( )
    {
        return compressionTime_ * 1.0e-2;
    }

    // Returns the downlink delay at the second receiving station in seconds.
    double getSecondReceivingStationDownlinkDelay( )
    {
        return secondReceivingStationDownlinkDelay_ * 1.0e-9;
    }

   private:
    int secondReceivingStationId_;
    int quasarOrSpacecraftId_;
    int phasePointIndicator_;
    int referenceFrequencyHighPart_;  // 2^24 mHz
    int referenceFrequencyLowPart_;   // mHz
    int composite1_;
    int compressionTime_;                      // 1e-2 sec
    int secondReceivingStationDownlinkDelay_;  // nsec
};

// TODO: test
// Derived class defining a Delta differential one-way ranging data block
class OdfDDorDataBlock : public OdfDataSpecificBlock
{
   public:
    /*!
     * Constructor. Parses a the observable specific portion of a Delta differential one-way ranging
     * block, according to table 3-4c of TRK-2-18 (2018).
     *
     * @param specificDataBits Observable specific portion of the data block.
     * @param dDorDataType Data type, specified according to section 3.2.4 of TRK-2-18 (2018). Can take values 5, 6.
     */
    OdfDDorDataBlock( const std::bitset< 128 > specificDataBits, const int dDorDataType );

    // Destructor
    ~OdfDDorDataBlock( ) { }

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

    // Returns the reference frequency in Hz
    double getReferenceFrequency( )
    {
        return std::pow( 2.0, 24 ) / 1.0E3 * referenceFrequencyHighPart_ +
                referenceFrequencyLowPart_ / 1.0E3;
    }

    int getComposite1( )
    {
        return composite1_;
    }

    // Returns the downlink delay at the second receiving station in seconds.
    double getSecondReceivingStationDownlinkDelay( )
    {
        return secondReceivingStationDownlinkDelay_ * 1.0e-9;
    }

   private:
    int secondReceivingStationId_;
    int quasarOrSpacecraftId_;
    int modulusIndicator_;
    int referenceFrequencyHighPart_;  // 2^24 mHz
    int referenceFrequencyLowPart_;   // mHz
    int composite1_;
    int modulusLowPart_;                       // 1e-7 sec
    int secondReceivingStationDownlinkDelay_;  // nsec
};

// Derived class defining an n-way (n >= 1) Doppler data block
class OdfDopplerDataBlock : public OdfDataSpecificBlock
{
   public:
    /*!
     * Constructor. Parses a the observable specific portion of an n-way Doppler block, according to
     * table 3-4d of TRK-2-18 (2018).
     *
     * @param specificDataBits Observable specific portion of the data block.
     * @param dopplerDataType Data type, specified according to section 3.2.4 of TRK-2-18 (2018). Can take values 11 (1-way),
     *      12 (2-way), 13 (3-way).
     */
    OdfDopplerDataBlock( const std::bitset< 128 > specificDataBits, const int dopplerDataType );

    // Destructor
    ~OdfDopplerDataBlock( ) { }

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

    // Returns the reference frequency in Hz
    double getReferenceFrequency( )
    {
        return std::pow( 2.0, 24 ) / 1.0E3 * referenceFrequencyHighPart_ +
                referenceFrequencyLowPart_ / 1.0E3;
    }

    // Returns the count interval in seconds.
    double getCompressionTime( )
    {
        return compressionTime_ * 1.0e-2;
    }

    // Returns the uplink delay at transmitting station in seconds.
    double getTransmittingStationUplinkDelay( )
    {
        return transmittingStationUplinkDelay_ * 1.0e-9;
    }

    /*!
     * Prints the data block to the specified output file.
     *
     * @param outFile File name on which to write data.
     */
    void printDataBlock( std::ofstream &outFile );

   private:
    int receiverChannel_;
    int spacecraftId_;
    int receiverExciterFlag_;
    int referenceFrequencyHighPart_;  // 2^24 mHz
    int referenceFrequencyLowPart_;   // mHz

    int reservedSegment_;
    int compressionTime_;  // 1e-2 sec

    int transmittingStationUplinkDelay_;  // nsec
};

// Derived class defining a sequential range data block
class OdfSequentialRangeDataBlock : public OdfDataSpecificBlock
{
   public:
    /*!
     * Constructor. Parses a the observable specific portion of an sequential range data block,
     * according to table 3-4e of TRK-2-18 (2018).
     *
     * @param dataBits Observable specific portion of the data block.
     */
    OdfSequentialRangeDataBlock( const std::bitset< 128 > dataBits );

    // Destructor
    ~OdfSequentialRangeDataBlock( ) { }

    int getSpacecraftId( )
    {
        return spacecraftId_;
    }

    // Returns the reference frequency in Hz
    double getReferenceFrequency( )
    {
        return std::pow( 2.0, 24 ) / 1.0E3 * referenceFrequencyHighPart_ +
                referenceFrequencyLowPart_ / 1.0E3;
    }

    // Returns the uplink delay at transmitting station in seconds.
    double getTransmittingStationUplinkDelay( )
    {
        return transmittingStationUplinkDelay_ * 1.0e-9;
    }

    // Returns the lowest ranging component
    int getLowestRangingComponent( )
    {
        return lowestRangingComponent_;
    }

    // Returns the uplink coder in-phase time offset in seconds
    int getUplinkCoderInPhaseTimeOffset( )
    {
        return uplinkCoderInPhaseTimeOffset_;
    }

    // Returns the composite two value in seconds
    double getCompositeTwo( )
    {
        return static_cast< double >( compositeTwo_ ) / 1.0e5;
    }

    int reservedBlock_;

   private:
    int spacecraftId_;

    int referenceFrequencyHighPart_;  // 2^24 mHz
    int referenceFrequencyLowPart_;   // mHz

    int transmittingStationUplinkDelay_;  // nsec

    int lowestRangingComponent_;
    int uplinkCoderInPhaseTimeOffset_;  // sec
    int compositeTwo_;                  // sec
};

// TODO: test
// Derived class defining a tone range data block
class OdfToneRangeDataBlock : public OdfDataSpecificBlock
{
   public:
    /*!
     * Constructor. Parses a the observable specific portion of a tone range data block, according
     * to table 3-4f of TRK-2-18 (2018).
     *
     * @param specificDataBits Observable specific portion of the data block.
     */
    OdfToneRangeDataBlock( const std::bitset< 128 > specificDataBits );

    // Destructor
    ~OdfToneRangeDataBlock( ) { }

    int getSpacecraftId( )
    {
        return spacecraftId_;
    }

    // Returns the reference frequency in Hz
    double getReferenceFrequency( )
    {
        return std::pow( 2.0, 24 ) / 1.0E3 * referenceFrequencyHighPart_ +
                referenceFrequencyLowPart_ / 1.0E3;
    }

    // Returns the uplink delay at transmitting station in seconds.
    double getTransmittingStationUplinkDelay( )
    {
        return transmittingStationUplinkDelay_ * 1.0e-9;
    }

   private:
    int integerObservableTime_;  // sec
    int spacecraftId_;
    int reservedBlock1_;

    int referenceFrequencyHighPart_;  // 2^24 mHz
    int referenceFrequencyLowPart_;   // mHz

    int reservedBlock2_;
    int reservedBlock3_;

    int transmittingStationUplinkDelay_;  // nsec
};

// TODO: test
// Derived class defining an angular observable data block
class OdfAngleDataBlock : public OdfDataSpecificBlock
{
   public:
    /*!
     * Constructor. Parses the observable specific portion of a angular data block, according to
     * table 3-4g of TRK-2-18 (2018).
     *
     * @param specificDataBits Observable specific portion of the data block.
     * @param angleDataType Data type, specified according to section 3.2.4 of TRK-2-18 (2018). Can take values in [51, 58].
     */
    OdfAngleDataBlock( const std::bitset< 128 > specificDataBits, const int angleDataType );

    // Destructor
    ~OdfAngleDataBlock( ) { }

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

// Class defining the common data data block of an ODF file, according to table 3-4a of TRK-2-18
// (2018).
class OdfCommonDataBlock
{
   public:
    /*!
     * Constructor. Parses the common portion of an ODF data block, according to table 3-4a of
     * TRK-2-18 (2018).
     *
     * @param commonDataBits Common portion of the data block.
     */
    OdfCommonDataBlock( const std::bitset< 160 > commonDataBits );

    // Returns the observable time in UTC seconds since the reference time specified in the header.
    Time getObservableTime( )
    {
        return Time( static_cast< double >( integerTimeTag_ ) ) +
                Time( static_cast< double >( fractionalTimeTag_ ) / 1000.0 );
    }

    // Returns the observable value in SI units.
    double getObservableValue( )
    {
        return static_cast< double >( integerObservable_ ) +
                static_cast< double >( fractionalObservable_ ) / 1.0E9;
    }

    // Returns the downlink delay at the receiving station in seconds.
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

    /*!
     * Prints the common data block to the specified output file.
     *
     * @param outFile File name on which to write data.
     */
    void printDataBlock( std::ofstream &outFile );

   private:
    uint32_t integerTimeTag_;            // sec
    int fractionalTimeTag_;              // msec
    int receivingStationDownlinkDelay_;  // nsec

    int integerObservable_;     // unit
    int fractionalObservable_;  // 1e-9 * unit
};

// Class defining a single data block of an ODF file.
class OdfDataBlock
{
   public:
    /*!
     * Constructor. Parses an ODF data block, dividing the data into a common data block and an
     * observable specific data block. According to section 3.2.4 of TRK-2-18 (2018).
     *
     * @param dataBits Data block of ODF file.
     */
    OdfDataBlock( const std::bitset< 288 > dataBits );

    // Returns the observable specific data block.
    std::shared_ptr< OdfDataSpecificBlock > getObservableSpecificDataBlock( )
    {
        return observableSpecificDataBlock_;
    }

    // Returns the common data block.
    std::shared_ptr< OdfCommonDataBlock > getCommonDataBlock( )
    {
        return commonDataBlock_;
    }

    /*!
     * Prints the data block to the specified output file.
     *
     * @param outFile File name on which to write data.
     */
    void printDataBlock( std::ofstream &outFile )
    {
        commonDataBlock_->printDataBlock( outFile );
        observableSpecificDataBlock_->printDataBlock( outFile );
        outFile << std::endl;
    }

   private:
    // Observable specific data block.
    std::shared_ptr< OdfDataSpecificBlock > observableSpecificDataBlock_;

    // Common data block.
    std::shared_ptr< OdfCommonDataBlock > commonDataBlock_;
};

// Class containing the raw data from an ODF file, according to TRK-2-18 (2018).
class OdfRawFileContents
{
   public:
    /*!
     * Constructor. Extracts all the data from an ODF file.
     *
     * @param odfFile File name/location of ODF file that is to be read
     */
    OdfRawFileContents( const std::string &odfFile );

    // File label group, table 3.2 of TRK-2-18 (2018)
    std::string systemId_;
    std::string programId_;
    uint32_t spacecraftId_;

    uint32_t fileCreationDate_;  // year, month, day (YYYMMDD): year from 1900
    uint32_t fileCreationTime_;  // hour, minute, second (HHMMSS)

    uint32_t fileReferenceDate_;  // year, month, day (YYYYMMDD)
    uint32_t fileReferenceTime_;  // hour, minute, second (HHMMSS)

    // ODF file name
    std::string fileName_;

    // Identifier group, table 3.3 of TRK-2-18 (2018)
    std::string identifierGroupStringA_;
    std::string identifierGroupStringB_;
    std::string identifierGroupStringC_;

    // Boolean indicating whether the EOF header was found (header should be present in all ODF
    // files)
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
    std::map< std::pair< int, int >, std::shared_ptr< OdfClockOffsetBlock > >
    getClockOffsetBlocks( )
    {
        return clockOffsetBlocks_;
    }

    /*!
     * Function to write the data in a binary ODF file to a txt file.
     *
     * @param odfTextFile Name of the output file.
     */
    void writeOdfToTextFile( const std::string &odfTextFile );

   private:
    //! Vector of data blocks
    std::vector< std::shared_ptr< OdfDataBlock > > dataBlocks_;

    //! Vector of ramp blocks indexed by transmitting station ID
    std::map< int, std::vector< std::shared_ptr< OdfRampBlock > > > rampBlocks_;

    //! Clock offset blocks indexed by pair of (primary station ID, secondary station ID)
    std::map< std::pair< int, int >, std::shared_ptr< OdfClockOffsetBlock > > clockOffsetBlocks_;

    /*!
     * Function to parse the contents of an ODF file label block, according to table 3.2 of TRK-2-18
     * (2018). The data contained in the block is returned by  reference.
     *
     * @param dataBits Data block (input).
     * @param systemId System ID (output).
     * @param programId Program ID (output)
     * @param spacecraftId Spacecraft ID (output)
     * @param fileCreationDate File creation date (output). Format: year, month, day (YYYMMDD): year from 1900
     * @param fileCreationTime File creation time (output). Format: hour, minute, second (HHMMSS)
     * @param fileReferenceDate File reference date (output). Format: year, month, day (YYYYMMDD)
     * @param fileReferenceTime File reference time (output). Format: hour, minute, second (HHMMSS)
     */
    void parseFileLabelData( std::bitset< 288 > dataBits,
                             std::string &systemId,
                             std::string &programId,
                             uint32_t &spacecraftId,
                             uint32_t &fileCreationDate,
                             uint32_t &fileCreationTime,
                             uint32_t &fileReferenceDate,
                             uint32_t &fileReferenceTime );

    /*!
     * Function to parse the contents of an ODF identifier block, according to table 3.3 of TRK-2-18
     * (2018). The data contained in the block is returned by  reference.
     *
     * @param dataBits Data block (input).
     * @param identifierGroupStringA identifierGroupStringA (output)
     * @param identifierGroupStringB identifierGroupStringB (output)
     * @param identifierGroupStringC identifierGroupStringC (output)
     */
    void parseIdentifierData( std::bitset< 288 > dataBits,
                              std::string &identifierGroupStringA,
                              std::string &identifierGroupStringB,
                              std::string &identifierGroupStringC );

    /*!
     * Function to parse the contents of an ODF file header block, according to table 3.1 of
     * TRK-2-18 (2018). The data contained in the header is returned by  reference. If the provided
     * dataBits block is not consistent with the format of a header, an error is thrown.
     *
     * @param dataBits Single header block of ODF file (input).
     * @param primaryKey Header primary key (output)
     * @param secondaryKey Header secondary key (output)
     * @param logicalRecordLength Logical record length (output)
     * @param groupStartPacketNumber Group start packet number (output)
     */
    void parseHeader( std::bitset< 288 > dataBits,
                      int32_t &primaryKey,
                      uint32_t &secondaryKey,
                      uint32_t &logicalRecordLength,
                      uint32_t &groupStartPacketNumber );

    /*!
     * Function to check if the current ODF data block is a header. Returns the primary key,
     * secondary key, logical record length, and group start packet number by reference if the group
     * IS a header.
     *
     * @param dataBits Single block of ODF file (input).
     * @param primaryKey Header primary key (output)
     * @param secondaryKey Header secondary key (output)
     * @param logicalRecordLength Logical record length (output)
     * @param groupStartPacketNumber Group start packet number (output)
     * @return Boolean indicating whether group is header or not.
     */
    bool currentBlockIsHeader( std::bitset< 288 > dataBits,
                               int &primaryKey,
                               unsigned int &secondaryKey,
                               unsigned int &logicalRecordLength,
                               unsigned int &groupStartPacketNumber );

    /*!
     * Function reads a single 36 byte block from ODF file. The block is returned by reference.
     *
     * @param file ODF file (input)
     * @param dataBits Read byte block (output).
     */
    void readOdfFileBlock( std::istream &file, std::bitset< 36 * 8 > &dataBits );
};

inline std::shared_ptr< OdfRawFileContents > readOdfFile( std::string fileName )
{
    return std::make_shared< OdfRawFileContents >( fileName );
}

}  // namespace input_output

}  // namespace tudat

#endif  // TUDAT_READ_ODF_FILE_H
