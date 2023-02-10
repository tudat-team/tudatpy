/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      The function printStandardScientificNotation() has been implemented to cope with
 *      cross-platform incompatibilities in the printed output of floating-point numbers in
 *      scientific notation.
 *
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

template< int FirstInputSize, int SecondInputSize >
std::bitset< FirstInputSize + SecondInputSize > mergeBitsets(
        std::bitset< FirstInputSize > firstInput,
        std::bitset< SecondInputSize > secondInput )
{
    std::bitset< FirstInputSize + SecondInputSize > returnBitset;
    for( int i = 0; i < FirstInputSize; i++ )
    {
        returnBitset[ i + SecondInputSize ] = firstInput[ i ];
    }

    for( int i = 0; i < SecondInputSize; i++ )
    {
        returnBitset[ i ] = secondInput[ i ];
    }
    return returnBitset;
}

template< int OutputBits, int InputBits >
std::bitset< OutputBits > getBitsetSegment(
        const std::bitset< InputBits > inputBits,
        const int startIndex )
{
    std::bitset< OutputBits > outputBits;

    // Check if final bit is valid
    if ( startIndex + OutputBits > InputBits )
    {
        throw std::runtime_error( "Error, when getting bit segment: requested bits are not part of the provided bitset." );
    }

    for( unsigned int i = 0; i < OutputBits; i++ )
    {
        outputBits[ i ] = inputBits[ InputBits - OutputBits - startIndex + i  ];
    }
    return outputBits;
}

template< int NumberOfBits >
int getSignedNBitInteger(
        std::bitset< NumberOfBits > inputBits )
{
    int outputInteger = -inputBits[ NumberOfBits - 1 ] * std::pow( 2, NumberOfBits - 1 );

    for( unsigned int i = 0; i < NumberOfBits - 1; i ++ )
    {
        outputInteger += inputBits[ i ] * std::pow( 2.0, i );
    }
    return outputInteger;
}

template < int NumberOfBytes >
void readBinaryFileBlock( std::istream& file,
                          std::bitset< NumberOfBytes * 8 >& dataBits )
{
    int numberOfBits = NumberOfBytes * 8;

    // TODO: Not being able to use read() with char dataChar [NumberOfBytes]... why?!
    char dataChar [NumberOfBytes][1];
    file.read( (char*)dataChar[0], NumberOfBytes );

    if ( !file.good( ) )
    {
        throw std::runtime_error( "Error when reading data block from ODF file." );
    }

    // Convert to bitset
    for ( int i = 0, bitCounter = 0; i < NumberOfBytes; ++i)
    {
        // Extract byte
        uint8_t byte = dataChar[0][i];
        for ( int j = 0; j < 8; ++j)
        {
            // Right shift byte to determine value of desired bit and save it to the bitset
            // Indexing of the byte and bitset starts from the right (i.e. the 0th bit is the rightmost one)
            dataBits[ numberOfBits - bitCounter - 1 ] = (byte >> ( 8 - 1 - j) ) & 1;
            ++bitCounter;
        }
    }
}

template< unsigned int NumberOfBits >
long convertBitsetToLong(const std::bitset< NumberOfBits >& bits) {
    if ( NumberOfBits > 32 )
    {
        throw std::runtime_error( "Error when converting bitset to long: specified number of bits (" +
            std::to_string(NumberOfBits) + "is larger than it is possible to represent with long (32).");
    }

    // Declare struct and create object s
    struct {
        // x with bit field of size numberOfBits
        // Sign extension is done automatically
        long x: NumberOfBits;
    } s;

    // Convert bitset to UNSIGNED long represented by numberOfBits bits. The remaining bits are determined by sign
    // extension, hence representing a SIGNED long.
    s.x = bits.to_ulong();
    return s.x;
}

template< unsigned int NumberBlockBits, unsigned int NumberItemBits, typename T >
void parseDataBlock (std::bitset< NumberBlockBits > dataBits,
                     const std::vector< bool >& unsignedItemFlag,
                     unsigned int argumentCounter,
                     unsigned int startBitCounter,
                     T& arg)
{
    if ( unsignedItemFlag.at( argumentCounter ) )
    {
        arg = getBitsetSegment< NumberItemBits, NumberBlockBits >( dataBits, startBitCounter ).to_ulong( );
    }
    else
    {
        arg = convertBitsetToLong< NumberItemBits >(
                getBitsetSegment< NumberItemBits, NumberBlockBits >( dataBits, startBitCounter ) );
    }

    ++argumentCounter;
    startBitCounter += NumberItemBits;

    if ( startBitCounter != NumberBlockBits )
    {
        throw std::runtime_error(
                "Error when parsing ODF file: block size (" + std::to_string( NumberBlockBits ) +
                " bits) and total item size (" + std::to_string(startBitCounter) + "bits ) are not consistent." );
    }
    else if ( argumentCounter != unsignedItemFlag.size() )
    {
        throw std::runtime_error(
                "Error when parsing ODF file: numbers of items (" + std::to_string( argumentCounter ) +
                ") and size of unsigned flag vector (" + std::to_string( unsignedItemFlag.size() ) +
                ") are not consistent." );
    }
}

template< unsigned int NumberBlockBits, unsigned int NumberItemBits, unsigned int... NumberItemBitsN,
        typename T, typename... TN >
void parseDataBlock (std::bitset< NumberBlockBits > dataBits,
                     const std::vector< bool >& unsignedItemFlag,
                     unsigned int argumentCounter,
                     unsigned int startBitCounter,
                     T& arg, TN&... args)
{
    if ( unsignedItemFlag.at( argumentCounter ) )
    {
        arg = getBitsetSegment< NumberItemBits, NumberBlockBits >( dataBits, startBitCounter ).to_ulong( );
    }
    else
    {
        arg = convertBitsetToLong< NumberItemBits >(
                getBitsetSegment< NumberItemBits, NumberBlockBits >( dataBits, startBitCounter ) );
    }

    ++argumentCounter;
    startBitCounter += NumberItemBits;

    parseDataBlock< NumberBlockBits, NumberItemBitsN ... >( dataBits, unsignedItemFlag, argumentCounter,
                                                            startBitCounter, args ...);
}

uint32_t convertCharactersToUnsignedInt32(
        unsigned char data[ 4 ] );

int32_t convertCharactersToSignedInt32(
        unsigned char data[ 4 ] );

class OdfClockOffsetBlock
{
public:

    unsigned int integerStartTime; // sec
    unsigned int fractionalStartTime; // nsec

    int integerClockOffset;  // sec
    int fractionalClockOffset; // nsec

    unsigned int primaryStationId;
    unsigned int secondaryStationId;

    unsigned int reservedBlock;

    unsigned int integerEndTime; // sec
    unsigned int fractionalEndTime; // nsec

    double getStartTime( )
    {
        return static_cast< double >( integerStartTime ) + static_cast< double >( fractionalStartTime ) * 1.0E-9;
    }

    double getEndTime( )
    {
        return static_cast< double >( integerEndTime ) + static_cast< double >( fractionalEndTime ) * 1.0E-9;
    }

    double getClockOffset( )
    {
        return static_cast< double >( integerClockOffset ) + static_cast< double >( fractionalClockOffset ) * 1.0E-9;
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
    std::map< int, std::vector< std::shared_ptr< OdfRampBlock > > > rampBlocks;
    std::shared_ptr< OdfClockOffsetBlock > clockOffsetBlock;
};

std::shared_ptr< OdfClockOffsetBlock > parseClockOffsetData( std::bitset< 288 > dataBits );

//! Function to parse the contents of an ODF orbit data block, specific for sequenctial range data.
std::shared_ptr< OdfSequentialRangeDataBlock > parseSequentialRangeData( unsigned char fileBlock[ 9 ][ 4 ], const int dopplerType );

std::shared_ptr< OdfSequentialRangeDataBlock > parseSequentialRangeData( std::bitset< 128 > dataBits );

//! Function to parse the contents of an ODF orbit data block, specific for Doppler data.
std::shared_ptr< OdfDopplerDataBlock > parseDopplerOrbitData( unsigned char fileBlock[ 9 ][ 4 ], const int dopplerType );

std::shared_ptr< OdfDopplerDataBlock > parseDopplerOrbitData( std::bitset< 128 > dataBits, const int dopplerType );

//! Function to parse the contents of an ODF orbit data block
std::shared_ptr< OdfDataBlock > parseOrbitData( unsigned char fileBlock[ 9 ][ 4 ] );

std::shared_ptr< OdfDataBlock > parseOrbitData( std::bitset< 288 > dataBits );

//! Function to parse the contents of an ODF ramp data block
OdfRampBlock parseRampData( unsigned char fileBlock[ 9 ][ 4 ] );

std::shared_ptr< OdfRampBlock > parseRampData( std::bitset< 288 > dataBits );

//! Function to parse the contents of an ODF file label block
void parseFileLabel( unsigned char fileBlock[ 9 ][ 4 ],
std::string& systemId, std::string& programId, std::string& fileCreationDate, std::string& fileCreationTme,
uint32_t& spacecraftIdNumber, uint32_t& fileReferenceDate, uint32_t& fileReferenceTime );

void parseFileLabelData(
        std::bitset< 288 > dataBits, std::string& systemId, std::string& programId, uint32_t& spacecraftId,
        uint32_t& fileCreationDate, uint32_t& fileCreationTime, uint32_t& fileReferenceDate,
        uint32_t& fileReferenceTime );

void parseIdentifierData(
        std::bitset< 288 > dataBits, std::string& identifierGroupStringA, std::string&identifierGroupStringB,
        std::string& identifierGroupStringC );

//! Function to parse the contents of an ODF file header block
void parseHeader( unsigned char fileBlock[ 9 ][ 4 ],
int32_t& primaryKey,
uint32_t& secondaryKey,
uint32_t& logicalrecordLength,
uint32_t& groupStartPacketNumber );

void parseHeader(
        std::bitset< 288 > dataBits, int32_t& primaryKey, uint32_t& secondaryKey, uint32_t& logicalRecordLength,
        uint32_t& groupStartPacketNumber);

//! Function to check if the current ODF data block is a header.
int currentBlockIsHeader(  unsigned char fileBlock[ 9 ][ 4 ], unsigned int& secondaryKeyInt );

bool currentBlockIsHeader(
        std::bitset< 288 > dataBits, int& primaryKey, unsigned int& secondaryKey,
        unsigned int& logicalRecordLength, unsigned int& groupStartPacketNumber );

//! Function to read a single 36 byte block from ODF file
void readOdfFileBlock(
        unsigned char fileBlock[ 9 ][ 4 ], std::istream& file );

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
