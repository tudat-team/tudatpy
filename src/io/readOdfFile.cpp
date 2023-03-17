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

    parseNumericalDataBlockWrapper< 288, 32, 32, 32, 32, 32, 32, 32, 32, 32 >(
            dataBits, unsignedCommonItemFlag,
            clockOffsetBlock->integerStartTime_,
            clockOffsetBlock->fractionalStartTime_,
            clockOffsetBlock->integerClockOffset_,
            clockOffsetBlock->fractionalClockOffset_,
            clockOffsetBlock->primaryStationId_,
            clockOffsetBlock->secondaryStationId_,
            clockOffsetBlock->reservedBlock_,
            clockOffsetBlock->integerEndTime_,
            clockOffsetBlock->fractionalEndTime_ );

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

    parseNumericalDataBlockWrapper< 128, 7, 10, 1, 22, 24, 20, 22, 22 >(
            dataBits, unsignedItemFlag,
            dopplerDataBlock->receiverChannel_,
            dopplerDataBlock->spacecraftId_,
            dopplerDataBlock->receiverExciterFlag_,
            dopplerDataBlock->referenceFrequencyHighPart_,
            dopplerDataBlock->referenceFrequencyLowPart_,
            dopplerDataBlock->reservedSegment_,
            dopplerDataBlock->compressionTime_,
            dopplerDataBlock->transmittingStationUplinkDelay_ );

//    std::cout << dopplerDataBlock->receiverChannel_ << " " <<
//              dopplerDataBlock->spacecraftId_ << " " <<
//              dopplerDataBlock->receiverExciterFlag_ << " " <<
//              std::fixed << dopplerDataBlock->getReferenceFrequency() << " " <<
//              dopplerDataBlock->reservedSegment_ << " " <<
//              dopplerDataBlock->compressionTime_ << " " <<
//              dopplerDataBlock->transmittingStationUplinkDelay_ << std::endl;

    return dopplerDataBlock;
}

// TODO: test
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

    parseNumericalDataBlockWrapper< 128, 7, 10, 1, 22, 24, 20, 22, 22 >(
            dataBits, unsignedItemFlag,
            rangeDataBlock->lowestRangingComponent_,
            rangeDataBlock->spacecraftId_,
            rangeDataBlock->reservedBlock_,
            rangeDataBlock->referenceFrequencyHighPart_,
            rangeDataBlock->referenceFrequencyLowPart_,
            rangeDataBlock->coderInPhaseTimeOffset_,
            rangeDataBlock->compositeTwo_,
            rangeDataBlock->transmittingStationUplinkDelay_ );

//    std::cout << rangeDataBlock->lowestRangingComponent_ << " " <<
//              rangeDataBlock->spacecraftId_ << " " <<
//              rangeDataBlock->reservedBlock_ << " " <<
//              std::fixed << rangeDataBlock->getReferenceFrequency() << " " <<
//              rangeDataBlock->coderInPhaseTimeOffset_ << " " <<
//              rangeDataBlock->compositeTwo_ << " " <<
//              rangeDataBlock->transmittingStationUplinkDelay_ << std::endl;

    return rangeDataBlock;
}

// TODO: test
std::shared_ptr< OdfDDodDataBlock > parseDDodOrbitData( std::bitset< 128 > dataBits, const int dDodType )
{
    std::shared_ptr< OdfDDodDataBlock > dataBlock = std::make_shared< OdfDDodDataBlock >( dDodType );

    // Items to read:
    // secondReceivingStationId_: uint7
    // quasarOrSpacecraftId_: uint10
    // phasePointIndicator_: uint1
    // referenceFrequencyHighPart_: uint22
    // referenceFrequencyLowPart_: uint24
    // composite1_: int20
    // compressionTime_: uint22
    // secondReceivingStationDownlinkDelay__: uint22

    std::vector< bool > unsignedItemFlag = { true, true, true, true, true, false, true, true };

    parseNumericalDataBlockWrapper< 128, 7, 10, 1, 22, 24, 20, 22, 22 >(
            dataBits, unsignedItemFlag,
            dataBlock->secondReceivingStationId_,
            dataBlock->quasarOrSpacecraftId_,
            dataBlock->phasePointIndicator_,
            dataBlock->referenceFrequencyHighPart_,
            dataBlock->referenceFrequencyLowPart_,
            dataBlock->composite1_,
            dataBlock->compressionTime_,
            dataBlock->secondReceivingStationDownlinkDelay_ );

    return dataBlock;
}

// TODO: test
std::shared_ptr< OdfDDorDataBlock > parseDDorOrbitData( std::bitset< 128 > dataBits, const int dDorType )
{
    std::shared_ptr< OdfDDorDataBlock > dataBlock = std::make_shared< OdfDDorDataBlock >( dDorType );

    // Items to read:
    // secondReceivingStationId_: uint7
    // quasarOrSpacecraftId_: uint10
    // modulusIndicator_: uint1
    // referenceFrequencyHighPart_: uint22
    // referenceFrequencyLowPart_: uint24
    // composite1_: int20
    // modulusLowPart_: uint22
    // secondReceivingStationDownlinkDelay__: uint22

    std::vector< bool > unsignedItemFlag = { true, true, true, true, true, false, true, true };

    parseNumericalDataBlockWrapper< 128, 7, 10, 1, 22, 24, 20, 22, 22 >(
            dataBits, unsignedItemFlag,
            dataBlock->secondReceivingStationId_,
            dataBlock->quasarOrSpacecraftId_,
            dataBlock->modulusIndicator_,
            dataBlock->referenceFrequencyHighPart_,
            dataBlock->referenceFrequencyLowPart_,
            dataBlock->composite1_,
            dataBlock->modulusLowPart_,
            dataBlock->secondReceivingStationDownlinkDelay_ );

    return dataBlock;
}

// TODO: test
std::shared_ptr< OdfToneRangeDataBlock > parseToneRangeOrbitData( std::bitset< 128 > dataBits )
{
    std::shared_ptr< OdfToneRangeDataBlock > dataBlock = std::make_shared< OdfToneRangeDataBlock >( );

    // Items to read:
    // integerObservableTime_: uint7
    // spacecraftId_: uint10
    // reservedBlock1_: uint1
    // referenceFrequencyHighPart_: uint22
    // referenceFrequencyLowPart_: uint24
    // reservedBlock2_: int20 -> incorrect information in TRK-2-18 2008
    // reservedBlock3_: uint22 -> incorrect information in TRK-2-18 2008
    // transmittingStationUplinkDelay_: uint22

    std::vector< bool > unsignedItemFlag = { true, true, true, true, true, false, true, true };

    parseNumericalDataBlockWrapper< 128, 7, 10, 1, 22, 24, 20, 22, 22 >(
            dataBits, unsignedItemFlag,
            dataBlock->integerObservableTime_,
            dataBlock->spacecraftId_,
            dataBlock->reservedBlock1_,
            dataBlock->referenceFrequencyHighPart_,
            dataBlock->referenceFrequencyLowPart_,
            dataBlock->reservedBlock2_,
            dataBlock->reservedBlock3_,
            dataBlock->transmittingStationUplinkDelay_ );

    return dataBlock;
}

// TODO: test
std::shared_ptr< OdfAngleDataBlock > parseAngleOrbitData( std::bitset< 128 > dataBits, const int angleType )
{
    std::shared_ptr< OdfAngleDataBlock > dataBlock = std::make_shared< OdfAngleDataBlock >( angleType );

    // Items to read:
    // reservedBlock1_: uint7
    // spacecraftId_: uint10
    // reservedBlock2_: uint1
    // reservedBlock3_: uint22
    // reservedBlock4_: uint24
    // reservedBlock5_: int20
    // reservedBlock6_: uint22
    // reservedBlock7_: uint22

    std::vector< bool > unsignedItemFlag = { true, true, true, true, true, false, true, true };

    parseNumericalDataBlockWrapper< 128, 7, 10, 1, 22, 24, 20, 22, 22 >(
            dataBits, unsignedItemFlag,
            dataBlock->reservedBlock1_,
            dataBlock->spacecraftId_,
            dataBlock->reservedBlock2_,
            dataBlock->reservedBlock3_,
            dataBlock->reservedBlock4_,
            dataBlock->reservedBlock5_,
            dataBlock->reservedBlock6_,
            dataBlock->reservedBlock7_ );

    return dataBlock;
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
    dataBlock->commonDataBlock_ = commonDataBlock;

    int dataType;

    std::bitset< 160 > commonDataBits = getBitsetSegment< 160, 288 >( dataBits, 0 );

    std::vector< bool > unsignedCommonItemFlag ( 14, true );
    unsignedCommonItemFlag.at( 3 ) = false;
    unsignedCommonItemFlag.at( 4 ) = false;

    parseNumericalDataBlockWrapper< 160, 32, 10, 22, 32, 32, 3, 7, 7, 2, 6, 2, 2, 2, 1 >(
            commonDataBits, unsignedCommonItemFlag,
            commonDataBlock->integerTimeTag_,
            commonDataBlock->fractionalTimeTag_,
            commonDataBlock->receivingStationDownlinkDelay_,
            commonDataBlock->integerObservable_,
            commonDataBlock->fractionalObservable_,
            commonDataBlock->formatId_,
            commonDataBlock->receivingStationId_,
            commonDataBlock->transmittingStationId_,
            commonDataBlock->transmittingStationNetworkId_,
            dataType,
            commonDataBlock->downlinkBandId_,
            commonDataBlock->uplinkBandId_,
            commonDataBlock->referenceBandId_,
            commonDataBlock->validity_ );

//    std::cout << commonDataBlock->integerTimeTag_ << " " <<
//        commonDataBlock->fractionalTimeTag_ << " " <<
//        commonDataBlock->receivingStationDownlinkDelay_ << " " <<
//                                                               commonDataBlock->integerObservable_ << " " <<
//                                                               commonDataBlock->fractionalObservable_ << " " <<
//                                                               commonDataBlock->formatId_ << " " <<
//                                                               commonDataBlock->receivingStationId_ << " " <<
//                                                               commonDataBlock->transmittingStationId_ << " " <<
//                                                               commonDataBlock->transmittingStationNetworkId_ << " " <<
//                                                               dataType << " " <<
//                                                               commonDataBlock->downlinkBandId_ << " " <<
//                                                               commonDataBlock->uplinkBandId_ << " " <<
//                                                               commonDataBlock->referenceBandId_ << " " <<
//                                                               commonDataBlock->validity_ << std::endl;

    // Read data type specific data
    std::bitset< 128 > specificDataBits = getBitsetSegment< 128, 288 >( dataBits, 160 );

    if ( dataType == 1 || dataType == 2 || dataType == 3 || dataType == 4 )
    {
        dataBlock->observableSpecificDataBlock_ = parseDDodOrbitData( specificDataBits, dataType );
    }
    else if ( dataType == 5 || dataType == 6 )
    {
        dataBlock->observableSpecificDataBlock_ = parseDDorOrbitData( specificDataBits, dataType );
    }
    else if ( dataType == 11 || dataType == 12 || dataType == 13 )
    {
        dataBlock->observableSpecificDataBlock_ = parseDopplerOrbitData( specificDataBits, dataType );
    }
    else if( dataType == 37 )
    {
        dataBlock->observableSpecificDataBlock_ = parseSequentialRangeData( specificDataBits );
    }
    else if ( dataType == 41 )
    {
        dataBlock->observableSpecificDataBlock_ = parseToneRangeOrbitData( specificDataBits );
    }
    else if ( dataType == 51 || dataType == 52 || dataType == 53 || dataType == 54 || dataType == 55 ||
                dataType == 56 || dataType == 57 || dataType == 58 )
    {
        dataBlock->observableSpecificDataBlock_ = parseAngleOrbitData( specificDataBits, dataType );
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

    parseNumericalDataBlockWrapper< 288, 32, 32, 32, 32, 22, 10, 32, 32, 32, 32 >(
            dataBits, unsignedCommonItemFlag,
            rampBlock->integerRampStartTime_,
            rampBlock->fractionalRampStartTime_,
            rampBlock->integerRampRate_,
            rampBlock->fractionalRampRate_,
            rampBlock->integerRampStartFrequency_,
            rampBlock->transmittingStationId_,
            rampBlock->integerRampStartFrequencyModulo_,
            rampBlock->fractionalRampStartFrequency_,
            rampBlock->integerRampEndTime_,
            rampBlock->fractionalRampEndTime_ );

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
    parseNumericalDataBlockWrapper< 4 * 32, 32, 32, 32, 32 >(
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
        parseNumericalDataBlockWrapper< 32, 32 >( testBits, { true }, testInteger );

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
    // systemId: ASCII-8
    // programId: ASCII-8
    const int numAsciiBits = 16 * 8, numAsciiBytes = 16;
    std::bitset< numAsciiBits > asciiDataBits = getBitsetSegment< numAsciiBits, 288 >( dataBits, 0 );

    parseStringsBlockWrapper< numAsciiBytes, 8, 8 >( asciiDataBits, systemId, programId );

    // Get bits with numerical data
    const int numNumericalBits = 20 * 8;
    std::bitset< numNumericalBits > numericalDataBits = getBitsetSegment< numNumericalBits, 288 >( dataBits, numAsciiBits );

    // Items to parse:
    // spacecraftId: uint32
    // fileCreationDate: uint32
    // fileCreationTime: uint32
    // fileReferenceDate: uint32
    // fileReferenceTime: uint32

    std::vector< bool > unsignedItemFlag = { true, true, true, true, true };

    // Read data from block
    parseNumericalDataBlockWrapper< numNumericalBits, 32, 32, 32, 32, 32 >(
            numericalDataBits, unsignedItemFlag,
            spacecraftId,
            fileCreationDate,
            fileCreationTime,
            fileReferenceDate,
            fileReferenceTime );
}

void parseIdentifierData( std::bitset< 288 > dataBits,
                          std::string& identifierGroupStringA,
                          std::string& identifierGroupStringB,
                          std::string& identifierGroupStringC )
{
    const int numBytes = 288 / 8;
    parseStringsBlockWrapper< numBytes, 8, 8, 20 >( dataBits, identifierGroupStringA, identifierGroupStringB,
                                                    identifierGroupStringC );
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
    if ( !dataFile.good( ) )
    {
        throw std::runtime_error( "Error when opening ODF file file." );
    }

    // Create OdfRawFileContents object, and add file name
    std::shared_ptr< OdfRawFileContents > odfFileContents = std::make_shared< OdfRawFileContents >( );
    odfFileContents->fileName_ = odfFile;

    // Declare used variables
    std::bitset< 36 * 8 > currenBitBlock;
    // Variables to parse headers
    int primaryKey;
    unsigned int secondaryKey, logicalRecordLength, groupStartPacketNumber;

    // Parse file label header
    readOdfFileBlock( dataFile, currenBitBlock );
    parseHeader( currenBitBlock, primaryKey, secondaryKey, logicalRecordLength, groupStartPacketNumber );

    if( primaryKey != 101 || secondaryKey != 0 || logicalRecordLength != 1 || groupStartPacketNumber != 0 )
    {
        throw std::runtime_error( "Error when reading ODF file, file label header invalid: primary key " +
        std::to_string( primaryKey ) + ", secondary key " + std::to_string( secondaryKey ) + ", logical record length" +
        std::to_string( logicalRecordLength ) + ", packet number " + std::to_string( groupStartPacketNumber ) + "." );
    }
//    std::cout << primaryKey << " " << secondaryKey << " " << logicalRecordLength << " " << groupStartPacketNumber << std::endl;

    // Parse file label data
    readOdfFileBlock( dataFile, currenBitBlock );
    parseFileLabelData( currenBitBlock, odfFileContents->systemId_, odfFileContents->programId_,
                        odfFileContents->spacecraftId_, odfFileContents->fileCreationDate_,
                        odfFileContents->fileCreationTime_,
                        odfFileContents->fileReferenceDate_, odfFileContents->fileReferenceTime_ );

    std::cout << odfFileContents->systemId_ << " " << odfFileContents->programId_ << " " << odfFileContents->spacecraftId_ <<
              " " << odfFileContents->fileCreationDate_ << " " << odfFileContents->fileCreationTime_ << " " <<
              odfFileContents->fileReferenceDate_ << " " << odfFileContents->fileReferenceTime_ << std::endl;

    // Parse identifier header
    readOdfFileBlock( dataFile, currenBitBlock );
    parseHeader( currenBitBlock, primaryKey, secondaryKey, logicalRecordLength, groupStartPacketNumber );

    if( primaryKey != 107 || secondaryKey != 0 || logicalRecordLength != 1 || groupStartPacketNumber != 2 )
    {
        throw std::runtime_error( "Error when reading ODF file, identifier header invalid: primary key " +
        std::to_string( primaryKey ) + ", secondary key " + std::to_string( secondaryKey ) + ", logical record length " +
        std::to_string( logicalRecordLength ) + ", packet number " + std::to_string( groupStartPacketNumber ) + "." );
    }

    // Parse identifier data
    readOdfFileBlock( dataFile, currenBitBlock );
    parseIdentifierData(
            currenBitBlock, odfFileContents->identifierGroupStringA_, odfFileContents->identifierGroupStringB_,
            odfFileContents->identifierGroupStringC_ );

    std::cout << odfFileContents->identifierGroupStringA_ << " " << odfFileContents->identifierGroupStringB_ << " " <<
              odfFileContents->identifierGroupStringC_ << std::endl;

    // Parse orbit data header
    readOdfFileBlock( dataFile, currenBitBlock );
    parseHeader( currenBitBlock, primaryKey, secondaryKey, logicalRecordLength, groupStartPacketNumber );
    if( primaryKey != 109 || secondaryKey != 0 || logicalRecordLength != 1 || groupStartPacketNumber != 4 )
    {
        throw std::runtime_error( "Error when reading ODF file, orbit header invalid: primary key " +
        std::to_string( primaryKey ) + ", secondary key " + std::to_string( secondaryKey ) + ", logical record length " +
        std::to_string( logicalRecordLength ) + ", packet number " + std::to_string( groupStartPacketNumber ) + "." );
    }

    // Reset data blocks and ramp blocks
    odfFileContents->dataBlocks_.resize( 0 );
    odfFileContents->rampBlocks_.clear();
    odfFileContents->clockOffsetBlocks_.clear();

    // Read file until summary or EOF header is found or eof() is reached
    for ( int currentRampStation = -1, currentBlockType = 109; !dataFile.eof( ); )
    {
        // Read current block
        readOdfFileBlock( dataFile, currenBitBlock );

        // Check if block is a header file
        bool blockIsHeader = currentBlockIsHeader(
                currenBitBlock, primaryKey, secondaryKey, logicalRecordLength, groupStartPacketNumber );

        // If block is header
        if ( blockIsHeader )
        {
            currentBlockType = primaryKey;
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
            }
            // Clock offset header
            else if ( primaryKey == 2040 )
            {
                if( primaryKey != 2040 || secondaryKey != 0 || logicalRecordLength != 1 )
                {
                    throw std::runtime_error( "Error when reading ODF file, clock offset header invalid: primary key " +
                    std::to_string( primaryKey ) + ", secondary key " + std::to_string( secondaryKey ) +
                    ", logical record length " + std::to_string( logicalRecordLength ) + "." );
                }
            }
            // Summary or EOF file header: exit loop
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
                std::shared_ptr< OdfDataBlock > currentDataBlock = parseOrbitData( currenBitBlock );
                if( currentDataBlock != nullptr )
                {
                    odfFileContents->dataBlocks_.push_back( currentDataBlock );
                }
            }
            // Read ramp data
            else if ( currentBlockType == 2030 )
            {
                std::shared_ptr< OdfRampBlock > currentDataBlock = parseRampData( currenBitBlock );
                if( currentDataBlock != nullptr )
                {
                    odfFileContents->rampBlocks_[ currentRampStation ].push_back( currentDataBlock );
                }
            }
            // Clock offset data
            else if ( currentBlockType == 2040 )
            {
                std::shared_ptr< OdfClockOffsetBlock > currentDataBlock = parseClockOffsetData( currenBitBlock );
                if( currentDataBlock != nullptr )
                {
                    odfFileContents->clockOffsetBlocks_[
                            std::make_pair( currentDataBlock->primaryStationId_,
                                            currentDataBlock->secondaryStationId_ ) ] = currentDataBlock;
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

    // Ignore summary data
    if ( primaryKey == 105 )
    {
        // Read summary data block
        readOdfFileBlock( dataFile, currenBitBlock );

        // Read next block
        readOdfFileBlock( dataFile, currenBitBlock );
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
        odfFileContents->eofHeaderFound_ = true;
    }
    else
    {
        odfFileContents->eofHeaderFound_ = false;
    }

    return odfFileContents;
}

} // namespace input_output

} // namespace tudat
