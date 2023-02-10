/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/readOdfFile.h"
#include "tudat/io/readBinaryFile.h"

namespace tudat
{
namespace input_output
{

// TODO: test
std::shared_ptr< OdfClockOffsetBlock > parseClockOffsetData( std::bitset< 288 > dataBits )
{

    std::shared_ptr< OdfClockOffsetBlock > clockOffsetBlock = std::make_shared< OdfClockOffsetBlock >( );

    // Data to parse
    // integerStartTime: uint32
    // fractionalStartTime: uint32
    // integerClockOffset: int32
    // fractionalClockOffset: int32
    // primaryStationId: uint32
    // secondaryStationId: uint32
    // reservedBlock: uint32
    // integerEndTime: uint32
    // fractionalEndTime: uint32

    std::vector< bool > unsignedCommonItemFlag ( 9, true );
    unsignedCommonItemFlag.at( 2 ) = false;
    unsignedCommonItemFlag.at( 3 ) = false;

    parseDataBlockWrapper< 288, 32, 32, 32, 32, 32, 32, 32, 32, 32 >(
            dataBits, unsignedCommonItemFlag,
            clockOffsetBlock->integerStartTime,
            clockOffsetBlock->fractionalStartTime,
            clockOffsetBlock->integerClockOffset,
            clockOffsetBlock->fractionalClockOffset,
            clockOffsetBlock->primaryStationId,
            clockOffsetBlock->secondaryStationId,
            clockOffsetBlock->reservedBlock,
            clockOffsetBlock->integerEndTime,
            clockOffsetBlock->fractionalEndTime );

    return clockOffsetBlock;
}

std::shared_ptr< OdfDopplerDataBlock > parseDopplerOrbitData( std::bitset< 128 > dataBits, const int dopplerType )
{
    std::shared_ptr< OdfDopplerDataBlock > dopplerDataBlock = std::make_shared< OdfDopplerDataBlock >( dopplerType );

    // Items to read:
    // receiverChannel: uint7
    // spacecraftId: uint10
    // receiverExciterFlag: uint1
    // referenceFrequencyHighPart: uint22
    // referenceFrequencyLowPart: uint24
    // reservedSegment: int20
    // compressionTime: uint22
    // transmittingStationDelay: uint22

    std::vector< bool > unsignedItemFlag = { true, true, true, true, true, false, true, true };

    parseDataBlockWrapper< 128, 7, 10, 1, 22, 24, 20, 22, 22 >(
            dataBits, unsignedItemFlag,
            dopplerDataBlock->receiverChannel,
            dopplerDataBlock->spacecraftId,
            dopplerDataBlock->receiverExciterFlag,
            dopplerDataBlock->referenceFrequencyHighPart,
            dopplerDataBlock->referenceFrequencyLowPart,
            dopplerDataBlock->reservedSegment,
            dopplerDataBlock->compressionTime,
            dopplerDataBlock->transmittingStationUplinkDelay );

    std::cout<< dopplerDataBlock->receiverChannel<<" "<<
                dopplerDataBlock->spacecraftId<<" "<<
                dopplerDataBlock->receiverExciterFlag<<" "<<
                std::fixed << dopplerDataBlock->getReferenceFrequency()<<" "<<
                dopplerDataBlock->reservedSegment<<" "<<
                dopplerDataBlock->compressionTime<<" "<<
                dopplerDataBlock->transmittingStationUplinkDelay<<std::endl;

    return dopplerDataBlock;
}

std::shared_ptr< OdfSequentialRangeDataBlock > parseSequentialRangeData( std::bitset< 128 > dataBits )
{
    std::shared_ptr< OdfSequentialRangeDataBlock > rangeDataBlock = std::make_shared< OdfSequentialRangeDataBlock >( );

    // Items to read:
    // lowestRangingComponent: uint7
    // spacecraftId: uint10
    // reservedBlock: uint1
    // referenceFrequencyHighPart: uint22
    // referenceFrequencyLowPart: uint24
    // coderInPhaseTimeOffset: int20
    // compositeTwo: uint22
    // transmittingStationUplinkDelay: uint22

    std::vector< bool > unsignedItemFlag = { true, true, true, true, true, false, true, true };

    parseDataBlockWrapper< 128, 7, 10, 1, 22, 24, 20, 22, 22 >(
            dataBits, unsignedItemFlag,
            rangeDataBlock->lowestRangingComponent,
            rangeDataBlock->spacecraftId,
            rangeDataBlock->reservedBlock,
            rangeDataBlock->referenceFrequencyHighPart,
            rangeDataBlock->referenceFrequencyLowPart,
            rangeDataBlock->coderInPhaseTimeOffset,
            rangeDataBlock->compositeTwo,
            rangeDataBlock->transmittingStationUplinkDelay );

    std::cout<< rangeDataBlock->lowestRangingComponent<<" "<<
            rangeDataBlock->spacecraftId<<" "<<
            rangeDataBlock->reservedBlock<<" "<<
            std::fixed << rangeDataBlock->getReferenceFrequency()<<" "<<
            rangeDataBlock->coderInPhaseTimeOffset<<" "<<
            rangeDataBlock->compositeTwo<<" "<<
            rangeDataBlock->transmittingStationUplinkDelay<<std::endl;

    return rangeDataBlock;
}

std::shared_ptr< OdfDataBlock > parseOrbitData( std::bitset< 288 > dataBits )
{
    std::shared_ptr< OdfDataBlock > dataBlock = std::make_shared< OdfDataBlock >( );

    // Read common data
    // integerTimeTag: uint32
    // fractionalTimeTag: uint10
    // receivingStationDownlinkDelay: uint22
    // integerObservable: int32
    // fractionalObservable: int32
    // formatId: uint3
    // receivingStation: uint7
    // transmittingStation: uint7
    // transmittingStationNetworkId: uint2
    // dataType: uint6
    // downlinkBand: uint2
    // uplinkBand: uint2
    // referenceBand: uint2
    // validity: uint1

    std::shared_ptr< OdfCommonDataBlock > commonDataBlock = std::make_shared< OdfCommonDataBlock >( );
    dataBlock->commonDataBlock = commonDataBlock;

    int dataType;

    std::bitset< 160 > commonDataBits = getBitsetSegment< 160, 288 >( dataBits, 0 );

    std::vector< bool > unsignedCommonItemFlag ( 14, true );
    unsignedCommonItemFlag.at( 3 ) = false;
    unsignedCommonItemFlag.at( 4 ) = false;

    parseDataBlockWrapper< 160, 32, 10, 22, 32, 32, 3, 7, 7, 2, 6, 2, 2, 2, 1 >(
            commonDataBits, unsignedCommonItemFlag,
            commonDataBlock->integerTimeTag,
            commonDataBlock->fractionalTimeTag,
            commonDataBlock->receivingStationDownlinkDelay,
            commonDataBlock->integerObservable,
            commonDataBlock->fractionalObservable,
            commonDataBlock->formatId,
            commonDataBlock->receivingStation,
            commonDataBlock->transmittingStation,
            commonDataBlock->transmittingStationNetworkId,
            dataType,
            commonDataBlock->downlinkBand,
            commonDataBlock->uplinkBand,
            commonDataBlock->referenceBand,
            commonDataBlock->validity );

            std::cout<< commonDataBlock->integerTimeTag<<" "<<
                        commonDataBlock->fractionalTimeTag<<" "<<
                        commonDataBlock->receivingStationDownlinkDelay<<" "<<
                        commonDataBlock->integerObservable<<" "<<
                        commonDataBlock->fractionalObservable<<" "<<
                        commonDataBlock->formatId<<" "<<
                        commonDataBlock->receivingStation<<" "<<
                        commonDataBlock->transmittingStation<<" "<<
                        commonDataBlock->transmittingStationNetworkId<<" "<<
                        dataType<<" "<<
                        commonDataBlock->downlinkBand<<" "<<
                        commonDataBlock->uplinkBand<<" "<<
                        commonDataBlock->referenceBand<<" "<<
                        commonDataBlock->validity<<std::endl;

    // Read data type specific data
    std::bitset< 128 > specificDataBits = getBitsetSegment< 128, 288 >( dataBits, 160 );

    if( dataType == 11 || dataType == 12 || dataType == 13 )
    {
        dataBlock->observableSpecificDataBlock = parseDopplerOrbitData( specificDataBits, dataType );
    }
    else if(  dataType == 37 )
    {
        dataBlock->observableSpecificDataBlock = parseSequentialRangeData( specificDataBits );
    }
    else
    {
        throw std::runtime_error( "Error, ODF data type " + std::to_string( dataType ) + " not recognized." );
    }

    return dataBlock;
}

std::shared_ptr< OdfRampBlock > parseRampData( std::bitset< 288 > dataBits )
{

    std::shared_ptr< OdfRampBlock > rampBlock = std::make_shared< OdfRampBlock >( );

    // Data to parse
    // integerRampStartTime: uint32
    // fractionalRampStartTime: uint32
    // integerRampRate: int32
    // fractionalRampRate: int32
    // integerRampStartFrequency: uint22
    // transmittingStationId: uint10
    // integerRampStartFrequencyModulo: uint32
    // fractionalRampStartFrequency: uint32
    // integerRampEndTime: uint32
    // fractionalRampEndTime: uint32

    std::vector< bool > unsignedCommonItemFlag ( 10, true );
    unsignedCommonItemFlag.at( 2 ) = false;
    unsignedCommonItemFlag.at( 3 ) = false;

    parseDataBlockWrapper< 288, 32, 32, 32, 32, 22, 10, 32, 32, 32, 32 >(
            dataBits, unsignedCommonItemFlag,
            rampBlock->integerRampStartTime,
            rampBlock->fractionalRampStartTime,
            rampBlock->integerRampRate,
            rampBlock->fractionalRampRate,
            rampBlock->integerRampStartFrequency,
            rampBlock->transmittingStationId,
            rampBlock->integerRampStartFrequencyModulo,
            rampBlock->fractionalRampStartFrequency,
            rampBlock->integerRampEndTime,
            rampBlock->fractionalRampEndTime );

//    std::cout<< rampBlock->transmittingStationId<<" "<<
//                std::fixed<<rampBlock->getRampStartTime()<<" "<<
//                std::fixed<<rampBlock->getRampRate()<<" "<<
//                std::fixed<<rampBlock->getRampStartFrequency()<<" "<<
//                std::fixed<<rampBlock->getRampEndTime()<<std::endl;
//    std::cout <<rampBlock->integerRampStartTime<<" "<<
//                rampBlock->fractionalRampStartTime<<std::endl;

    return rampBlock;
}

void parseHeader( std::bitset< 288 > dataBits,
                  int32_t& primaryKey,
                  uint32_t& secondaryKey,
                  uint32_t& logicalRecordLength,
                  uint32_t& groupStartPacketNumber)
{
    // Items to read:
    // primaryKey: int32
    // secondaryKey: uint32
    // logicalRecordLength: uint32
    // groupStartPacketNumber: uint32
    // filler: uint32 * 5

    std::vector< bool > unsignedItemFlag = { false, true, true, true };
    std::bitset< 4 * 32 > validDataBits = getBitsetSegment< 4 * 32, 288 >( dataBits, 0 );

    // Read data from block
    parseDataBlockWrapper< 4 * 32, 32, 32, 32, 32 >(
            validDataBits, unsignedItemFlag,
            primaryKey,
            secondaryKey,
            logicalRecordLength,
            groupStartPacketNumber );

    // Read subsequent blocks and check if zero
    uint32_t testInteger;
    for( int i = 0; i < 5; i++ )
    {
        std::bitset< 32 > testBits = getBitsetSegment< 32, 288 >( dataBits, ( i + 4 ) * 32 );
        parseDataBlockWrapper< 32, 32 >( testBits, { true },  testInteger );

        if( testInteger != 0 )
        {
            throw std::runtime_error( "Error when reading ODF file, header file inconsistent: filler items are not zero." );
        }
    }
}

void parseFileLabelData(std::bitset< 288 > dataBits,
                        std::string& systemId,
                        std::string& programId,
                        uint32_t& spacecraftId,
                        uint32_t& fileCreationDate,
                        uint32_t& fileCreationTime,
                        uint32_t& fileReferenceDate,
                        uint32_t& fileReferenceTime )
{
    // Get bits with ascii data
    std::bitset< 16 * 8 > asciiDataBits = getBitsetSegment< 16 * 8, 288 >( dataBits, 0 );

    // systemId: ASCII-8
    systemId.resize( 8 );
    // programId: ASCII-8
    programId.resize( 8 );

    for ( int i = 0; i < 8; ++i )
    {
        systemId[i] = getBitsetSegment< 8, 16 * 8 >( asciiDataBits, i * 8 ).to_ulong();
        programId[i] = getBitsetSegment< 8, 16 * 8 >( asciiDataBits, 8 * 8 + i * 8 ).to_ulong();
    }

    // Get bits with numerical data
    std::bitset< 20 * 8 > numericalDataBits = getBitsetSegment< 20 * 8, 288 >( dataBits, 16 * 8 );

    // Items to parse:
    // spacecraftId: uint32
    // fileCreationDate: uint32
    // fileCreationTime: uint32
    // fileReferenceDate: uint32
    // fileReferenceTime: uint32

    std::vector< bool > unsignedItemFlag = { true, true, true, true, true };

    // Read data from block
    parseDataBlockWrapper< 5 * 32, 32, 32, 32, 32, 32 >(
            numericalDataBits, unsignedItemFlag,
            spacecraftId,
            fileCreationDate,
            fileCreationTime,
            fileReferenceDate,
            fileReferenceTime );
}

void parseIdentifierData( std::bitset< 288 > dataBits,
                          std::string& identifierGroupStringA,
                          std::string&identifierGroupStringB,
                          std::string& identifierGroupStringC )
{
    // Read three strings from data blocks
    identifierGroupStringA.resize( 8 );
    identifierGroupStringB.resize( 8 );
    identifierGroupStringC.resize( 20 );

    for ( int i = 0; i < 8; ++i )
    {
        identifierGroupStringA[i] = getBitsetSegment< 8, 288 >( dataBits, i * 8 ).to_ulong();
        identifierGroupStringB[i] = getBitsetSegment< 8, 288 >( dataBits, 8 * 8 + i * 8 ).to_ulong();
    }

    for ( int i = 0; i < 20; ++i )
    {
        identifierGroupStringC[i] = getBitsetSegment< 8, 288 >( dataBits, 16 * 8 + i * 8 ).to_ulong();
    }
}

bool currentBlockIsHeader( std::bitset< 288 > dataBits,
                          int& primaryKey,
                          unsigned int& secondaryKey,
                          unsigned int& logicalRecordLength,
                          unsigned int& groupStartPacketNumber )
{
    bool blockIsHeader;

    // Try reading current block as a header, if an exception is caught, block is not a header.
    try
    {
        parseHeader( dataBits, primaryKey, secondaryKey, logicalRecordLength,
                     groupStartPacketNumber );

        blockIsHeader =  true;
    }
    catch( std::runtime_error )
    {
        blockIsHeader =  false;
    }
    return blockIsHeader;
}

void readOdfFileBlock( std::istream& file,
                       std::bitset< 36 * 8 >& dataBits )
{
    readBinaryFileBlock< 36 >( file, dataBits );
}

//! Function to read the contents of an ODF file into an OdfRawFileContents object
std::shared_ptr< OdfRawFileContents > readOdfFile(
        const std::string& odfFile )
{
    // Open file
    std::ifstream dataFile( odfFile, std::ios_base::binary);

    // Create OdfRawFileContents object, and add file name
    std::shared_ptr< OdfRawFileContents > odfFileContents = std::make_shared< OdfRawFileContents >( );
    odfFileContents->fileName = odfFile;

    // Declare used variables
    std::bitset< 36 * 8 > currentFileBlock;
    // Variables to parse headers
    int primaryKey;
    unsigned int secondaryKey, logicalRecordLength, groupStartPacketNumber;

    // Parse file label header
    readOdfFileBlock( dataFile, currentFileBlock );
    parseHeader( currentFileBlock, primaryKey, secondaryKey, logicalRecordLength, groupStartPacketNumber );

    if( primaryKey != 101 || secondaryKey != 0 || logicalRecordLength != 1 || groupStartPacketNumber != 0 )
    {
        throw std::runtime_error( "Error when reading ODF file, file label header invalid: primary key " +
        std::to_string( primaryKey ) + ", secondary key " + std::to_string( secondaryKey ) + ", logical record length" +
        std::to_string( logicalRecordLength ) + ", packet number " + std::to_string( groupStartPacketNumber ) + "." );
    }
//    std::cout << primaryKey << " " << secondaryKey << " " << logicalRecordLength << " " << groupStartPacketNumber << std::endl;

    // Parse file label data
    readOdfFileBlock( dataFile, currentFileBlock );
    parseFileLabelData( currentFileBlock, odfFileContents->systemId, odfFileContents->programId,
                        odfFileContents->spacecraftId, odfFileContents->fileCreationDate,
                        odfFileContents->fileCreationTime,
                        odfFileContents->fileReferenceDate, odfFileContents->fileReferenceTime );

    std::cout<< odfFileContents->systemId << " " << odfFileContents->programId << " " << odfFileContents->spacecraftId <<
        " " << odfFileContents->fileCreationDate << " " << odfFileContents->fileCreationTime << " " <<
        odfFileContents->fileReferenceDate << " " << odfFileContents->fileReferenceTime << std::endl;

    // Parse identifier header
    readOdfFileBlock( dataFile, currentFileBlock );
    parseHeader( currentFileBlock, primaryKey, secondaryKey, logicalRecordLength, groupStartPacketNumber );

    if( primaryKey != 107 || secondaryKey != 0 || logicalRecordLength != 1 || groupStartPacketNumber != 2 )
    {
        throw std::runtime_error( "Error when reading ODF file, identifier header invalid: primary key " +
        std::to_string( primaryKey ) + ", secondary key " + std::to_string( secondaryKey ) + ", logical record length " +
        std::to_string( logicalRecordLength ) + ", packet number " + std::to_string( groupStartPacketNumber ) + "." );
    }

    // Parse identifier data
    readOdfFileBlock( dataFile, currentFileBlock );
    parseIdentifierData(
            currentFileBlock, odfFileContents->identifierGroupStringA, odfFileContents->identifierGroupStringB,
            odfFileContents->identifierGroupStringC );

    std::cout<< odfFileContents->identifierGroupStringA << " " << odfFileContents->identifierGroupStringB << " " <<
        odfFileContents->identifierGroupStringC << std::endl;

    // Parse orbit data header
    readOdfFileBlock( dataFile, currentFileBlock );
    parseHeader( currentFileBlock, primaryKey, secondaryKey, logicalRecordLength, groupStartPacketNumber );
    if( primaryKey != 109 || secondaryKey != 0 || logicalRecordLength != 1 || groupStartPacketNumber != 4 )
    {
        throw std::runtime_error( "Error when reading ODF file, orbit header invalid: primary key " +
        std::to_string( primaryKey ) + ", secondary key " + std::to_string( secondaryKey ) + ", logical record length " +
        std::to_string( logicalRecordLength ) + ", packet number " + std::to_string( groupStartPacketNumber ) + "." );
    }

    // Reset data blocks and ramp blocks
    odfFileContents->dataBlocks.resize( 0 );
    odfFileContents->rampBlocks.clear();

    // Read file until a non-ramp header is found or EOF is reached
    for ( int currentRampStation = -1, currentBlockType = 109; !dataFile.eof( ); )
    {
        // Read current block
        readOdfFileBlock( dataFile, currentFileBlock );

        // Check if block is a header file
        bool blockIsHeader = currentBlockIsHeader(
                currentFileBlock, primaryKey, secondaryKey, logicalRecordLength, groupStartPacketNumber );

        // If block is header
        if ( blockIsHeader )
        {
            // Ramp group header
            if ( primaryKey == 2030 )
            {
//                std::cout<<primaryKey<<" "<<secondaryKey<<" "<<logicalRecordLength<<" "<<groupStartPacketNumber<<std::endl;
                if( secondaryKey < 0 || secondaryKey > 99 || logicalRecordLength != 1 )
                {
                    throw std::runtime_error( "Error when reading ODF file, ramp header invalid: primary key " +
                    std::to_string( primaryKey ) + ", secondary key " + std::to_string( secondaryKey ) +
                    ", logical record length " + std::to_string( logicalRecordLength ) + "." );
                }
                currentRampStation = secondaryKey;
                currentBlockType = primaryKey;
            }
            // Clock offset, summary, or EOF file header: exit loop
            else
            {
                break;
            }
        }
        // If not a  header, read associated block
        else
        {
            // Read orbit data
            if ( currentBlockType == 109 )
            {
                std::shared_ptr< OdfDataBlock > currentDataBlock = parseOrbitData( currentFileBlock );
                if( currentDataBlock != nullptr )
                {
                    odfFileContents->dataBlocks.push_back( currentDataBlock );
                }
            }
            // Read ramp data
            else if ( currentBlockType == 2030 )
            {
                std::shared_ptr< OdfRampBlock > currentDataBlock = parseRampData( currentFileBlock );
                if( currentDataBlock != nullptr )
                {
                    odfFileContents->rampBlocks[ currentRampStation ].push_back( currentDataBlock );
                }
            }
            else
            {
                throw std::runtime_error( "Error when reading ODF group, invalid block type." );
            }
        }
    }

    if ( dataFile.eof( ) )
    {
        throw std::runtime_error( "Error when reading ODF file: end of file was found before EOF group." );
    }

    // Parse clock offset data
    if ( primaryKey == 2040 )
    {
        if( primaryKey != 2040 || secondaryKey != 0 || logicalRecordLength != 1 )
        {
            throw std::runtime_error( "Error when reading ODF file, clock offset header invalid: primary key " +
            std::to_string( primaryKey ) + ", secondary key " + std::to_string( secondaryKey ) + ", logical record length " +
            std::to_string( logicalRecordLength ) + "." );
        }

        readOdfFileBlock( dataFile, currentFileBlock );
        odfFileContents->clockOffsetBlock = parseClockOffsetData( currentFileBlock );

        // Read next block
        readOdfFileBlock( dataFile, currentFileBlock );
    }

    // Ignore summary data
    if ( primaryKey == 105 )
    {
        // Read summary data block
        readOdfFileBlock( dataFile, currentFileBlock );

        // Read next block
        readOdfFileBlock( dataFile, currentFileBlock );
    }

    // EOF group
    if ( primaryKey == -1 )
    {
        if( primaryKey != -1 || secondaryKey != 0 || logicalRecordLength != 0 )
        {
            throw std::runtime_error( "Error when reading ODF file, EOF header invalid: primary key " +
            std::to_string( primaryKey ) + ", secondary key " + std::to_string( secondaryKey ) + ", logical record length " +
            std::to_string( logicalRecordLength ) + "." );
        }
        odfFileContents->eofHeaderFound = true;
    }
    else
    {
        odfFileContents->eofHeaderFound = false;
    }

    return odfFileContents;
}

} // namespace input_output

} // namespace tudat
