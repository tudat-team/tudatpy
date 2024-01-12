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

namespace tudat
{
namespace input_output
{

OdfClockOffsetBlock::OdfClockOffsetBlock( const std::bitset< 288 > dataBits )
{
    // Data to parse:
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
            integerStartTime_,
            fractionalStartTime_,
            integerClockOffset_,
            fractionalClockOffset_,
            primaryStationId_,
            secondaryStationId_,
            reservedBlock_,
            integerEndTime_,
            fractionalEndTime_ );
}

void OdfRampBlock::printDataBlock( std::ofstream& outFile )
{
    outFile << std::setfill(' ') << std::setw(3) << transmittingStationId_ << " ";
    outFile << std::setfill(' ') << std::setw(20) << std::setprecision(9) << getRampStartTime() << " ";
    outFile << std::setfill(' ') << std::setw(20) << std::setprecision(9) << getRampRate() << " ";
    outFile << std::setfill(' ') << std::setw(20) << std::setprecision(9) << getRampStartFrequency() << " ";
    outFile << std::setfill(' ') << std::setw(20) << std::setprecision(9) << getRampEndTime() << " ";
    outFile << std::endl;
}


OdfDDodDataBlock::OdfDDodDataBlock( const std::bitset< 128 > specificDataBits, const int dDodDataType ):
    OdfDataSpecificBlock( dDodDataType )
{
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
            specificDataBits, unsignedItemFlag,
            secondReceivingStationId_,
            quasarOrSpacecraftId_,
            phasePointIndicator_,
            referenceFrequencyHighPart_,
            referenceFrequencyLowPart_,
            composite1_,
            compressionTime_,
            secondReceivingStationDownlinkDelay_ );
}

OdfDDorDataBlock::OdfDDorDataBlock( const std::bitset< 128 > specificDataBits, const int dDorDataType ):
    OdfDataSpecificBlock( dDorDataType )
{
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
            specificDataBits, unsignedItemFlag,
            secondReceivingStationId_,
            quasarOrSpacecraftId_,
            modulusIndicator_,
            referenceFrequencyHighPart_,
            referenceFrequencyLowPart_,
            composite1_,
            modulusLowPart_,
            secondReceivingStationDownlinkDelay_ );
}

OdfDopplerDataBlock::OdfDopplerDataBlock(
        const std::bitset< 128 > specificDataBits,
        const int dopplerDataType ):
    OdfDataSpecificBlock( dopplerDataType )
{
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
            specificDataBits, unsignedItemFlag,
            receiverChannel_,
            spacecraftId_,
            receiverExciterFlag_,
            referenceFrequencyHighPart_,
            referenceFrequencyLowPart_,
            reservedSegment_,
            compressionTime_,
            transmittingStationUplinkDelay_ );
}

void OdfDopplerDataBlock::printDataBlock( std::ofstream& outFile )
{
    outFile << std::setfill(' ') << std::setw(2) << receiverChannel_ << " ";
    outFile << std::setfill(' ') << std::setw(3) << spacecraftId_ << " ";
    outFile << std::setfill(' ') << std::setw(2) << receiverExciterFlag_ << " ";
    outFile << std::setfill(' ') << std::setw(14) << std::setprecision(3) << getReferenceFrequency() << " ";
    outFile << std::setfill(' ') << std::setw(7) << reservedSegment_ << " ";
    outFile << std::setfill(' ') << std::setw(7) << std::setprecision(3) << getCompressionTime() << " ";
    outFile << std::setfill(' ') << std::setw(7) << std::setprecision(3) << transmittingStationUplinkDelay_ << " ";
}

OdfSequentialRangeDataBlock::OdfSequentialRangeDataBlock( const std::bitset< 128 > dataBits ):
    OdfDataSpecificBlock( 37 )
{
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
            lowestRangingComponent_,
            spacecraftId_,
            reservedBlock_,
            referenceFrequencyHighPart_,
            referenceFrequencyLowPart_,
            uplinkCoderInPhaseTimeOffset_,
            compositeTwo_,
            transmittingStationUplinkDelay_ );
}

OdfToneRangeDataBlock::OdfToneRangeDataBlock( const std::bitset< 128 > specificDataBits ):
    OdfDataSpecificBlock( 41 )
{
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
            specificDataBits, unsignedItemFlag,
            integerObservableTime_,
            spacecraftId_,
            reservedBlock1_,
            referenceFrequencyHighPart_,
            referenceFrequencyLowPart_,
            reservedBlock2_,
            reservedBlock3_,
            transmittingStationUplinkDelay_ );
}

OdfAngleDataBlock::OdfAngleDataBlock( const std::bitset< 128 > specificDataBits, const int angleDataType ):
    OdfDataSpecificBlock( angleDataType )
{
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
            specificDataBits, unsignedItemFlag,
            reservedBlock1_,
            spacecraftId_,
            reservedBlock2_,
            reservedBlock3_,
            reservedBlock4_,
            reservedBlock5_,
            reservedBlock6_,
            reservedBlock7_ );
}

OdfCommonDataBlock::OdfCommonDataBlock( const std::bitset< 160 > commonDataBits )
{
    // Items to read:
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

    std::vector< bool > unsignedCommonItemFlag ( 14, true );
    unsignedCommonItemFlag.at( 3 ) = false;
    unsignedCommonItemFlag.at( 4 ) = false;

    parseNumericalDataBlockWrapper< 160, 32, 10, 22, 32, 32, 3, 7, 7, 2, 6, 2, 2, 2, 1 >(
            commonDataBits, unsignedCommonItemFlag,
            integerTimeTag_,
            fractionalTimeTag_,
            receivingStationDownlinkDelay_,
            integerObservable_,
            fractionalObservable_,
            formatId_,
            receivingStationId_,
            transmittingStationId_,
            transmittingStationNetworkId_,
            dataType_,
            downlinkBandId_,
            uplinkBandId_,
            referenceBandId_,
            validity_ );

    if ( formatId_ != 2 )
    {
        throw std::runtime_error( "Error when reading ODF file: reading of ODF files with format ID " + std::to_string( formatId_ ) +
            " not implemented." );
    }
}

void OdfCommonDataBlock::printDataBlock( std::ofstream& outFile )
{
    outFile << std::fixed << std::setprecision(3) << getObservableTime() << " ";
    outFile << std::setfill(' ') << std::setw(11) << receivingStationDownlinkDelay_ << " ";
    outFile << std::setfill(' ') << std::setw(20) << std::setprecision(9) << getObservableValue() << " ";
    outFile << std::setfill(' ') << std::setw(2) << formatId_ << " ";
    outFile << std::setfill(' ') << std::setw(4) << receivingStationId_ << " ";
    outFile << std::setfill(' ') << std::setw(4) << transmittingStationId_ << " ";
    outFile << std::setfill(' ') << std::setw(3) << transmittingStationNetworkId_ << " ";
    outFile << std::setfill(' ') << std::setw(5) << dataType_ << " ";
    outFile << std::setfill(' ') << std::setw(3) << downlinkBandId_ << " ";
    outFile << std::setfill(' ') << std::setw(2) << uplinkBandId_ << " ";
    outFile << std::setfill(' ') << std::setw(2) << referenceBandId_ << " ";
    outFile << std::setfill(' ') << std::setw(1) << validity_ << " ";
}

OdfDataBlock::OdfDataBlock( std::bitset< 288 > dataBits )
{
    // Extract bitset with common data
    std::bitset< 160 > commonDataBits = getBitsetSegment< 160, 288 >( dataBits, 0 );
    // Extract bitset with data-type specific data
    std::bitset< 128 > specificDataBits = getBitsetSegment< 128, 288 >( dataBits, 160 );

    commonDataBlock_ = std::make_shared< OdfCommonDataBlock >( commonDataBits );

    int dataType = commonDataBlock_->dataType_;
    if ( dataType == 1 || dataType == 2 || dataType == 3 || dataType == 4 )
    {
        observableSpecificDataBlock_ = std::make_shared< OdfDDodDataBlock >( specificDataBits, dataType );
    }
    else if ( dataType == 5 || dataType == 6 )
    {
        observableSpecificDataBlock_ = std::make_shared< OdfDDorDataBlock >( specificDataBits, dataType );
    }
    else if ( dataType == 11 || dataType == 12 || dataType == 13 )
    {
        observableSpecificDataBlock_ = std::make_shared< OdfDopplerDataBlock >( specificDataBits, dataType );
    }
    else if( dataType == 37 )
    {
        observableSpecificDataBlock_ = std::make_shared< OdfSequentialRangeDataBlock >( specificDataBits );
    }
    else if ( dataType == 41 )
    {
        observableSpecificDataBlock_ = std::make_shared< OdfToneRangeDataBlock >( specificDataBits );
    }
    else if ( dataType == 51 || dataType == 52 || dataType == 53 || dataType == 54 || dataType == 55 ||
                dataType == 56 || dataType == 57 || dataType == 58 )
    {
        observableSpecificDataBlock_ = std::make_shared< OdfAngleDataBlock >( specificDataBits, dataType );
    }
    else
    {
        throw std::runtime_error( "Error, ODF data type " + std::to_string( dataType ) + " not recognized." );
    }
}

OdfRampBlock::OdfRampBlock( const std::bitset< 288 > dataBits )
{
    // Data to parse:
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
            integerRampStartTime_,
            fractionalRampStartTime_,
            integerRampRate_,
            fractionalRampRate_,
            integerRampStartFrequency_,
            transmittingStationId_,
            integerRampStartFrequencyModulo_,
            fractionalRampStartFrequency_,
            integerRampEndTime_,
            fractionalRampEndTime_ );
}

void OdfRawFileContents::parseHeader( std::bitset< 288 > dataBits,
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

void OdfRawFileContents::writeOdfToTextFile( const std::string& odfTextFile )
{
    std::ofstream dataFile( odfTextFile );
    if ( !dataFile.good( ) )
    {
        throw std::runtime_error( "Error when opening output ODF text file." );
    }

    dataFile << "================================================" << std::endl;
    dataFile << "---Label Group:        Header Record------------" << std::endl;
    dataFile << "      Prime_Key    Second_ Key    Log_Rec_Len   Gp_St_Pkt_No" << std::endl;
    dataFile << "            101              0              1              0" << std::endl;
    dataFile << "-----------------------Data---------------------" << std::endl;
    dataFile << "  SystemID  ProgrmID  SCID  CrDate  CrTime  FRefDate  FRefTime" << std::endl;
    dataFile << "  " << systemId_ << "  " << programId_ << "   " << spacecraftId_ << "  " << fileCreationDate_ << "  " << fileCreationTime_ <<
        "  " << fileReferenceDate_ << "    " << std::setw(6) << std::setfill('0') << fileReferenceTime_ << std::endl;

    dataFile << "================================================" << std::endl;
    dataFile << "---Identfier Group:    Header Record------------" << std::endl;
    dataFile << "      Prime_Key    Second_ Key    Log_Rec_Len   Gp_St_Pkt_No" << std::endl;
    dataFile << "            107              0              1              2" << std::endl;

    dataFile << "-----------------------Data---------------------" << std::endl;
    dataFile << "   Item_1    Item_2          Item_3" << std::endl;
    dataFile << "  " << identifierGroupStringA_ << "  " << identifierGroupStringB_ << "   " << identifierGroupStringC_ << std::endl;

    dataFile << "================================================" << std::endl;
    dataFile << "---Orbit Data Group:   Header Record------------" << std::endl;
    dataFile << "      Prime_Key    Second_ Key    Log_Rec_Len   Gp_St_Pkt_No" << std::endl;
    dataFile << "            109              0              1              4" << std::endl;

    dataFile << "-----------------------Data---------------------" << std::endl;
    dataFile << "   Time Tag     dl_delay          Observable   Fmt  DSSr DSSt Net D_Typ DL UL Ex V 15  16 17 Reference Freq Item_20 Item_21 Item_22" << std::endl;


    int packetCounter = 4;
    for ( unsigned int i = 0; i < dataBlocks_.size(); ++i )
    {
        ++packetCounter;
        dataBlocks_.at( i )->printDataBlock( dataFile );
    }

    for ( auto it = rampBlocks_.begin(); it != rampBlocks_.end(); ++it )
    {
        ++packetCounter;
        dataFile << "================================================" << std::endl;
        dataFile << "---Ramp Group:         Header Record------------" << std::endl;
        dataFile << "      Prime_Key    Second_ Key    Log_Rec_Len   Gp_St_Pkt_No" << std::endl;
        dataFile << "             -1    " << std::setw(11) << it->first << "              1   "
            << std::setw(12) << packetCounter << std::endl;
        dataFile << "-----------------------Data---------------------" << std::endl;
        dataFile << "DSS   Ramp Start Time         Ramp Rate         Start Frequency        Ramp End Time" << std::endl;

        for ( unsigned int i = 0; i < it->second.size(); ++i )
        {
            ++packetCounter;
            it->second.at( i )->printDataBlock( dataFile );
        }
    }

    if ( eofHeaderFound_ )
    {
        ++packetCounter;
        dataFile << "================================================" << std::endl;
        dataFile << "      Prime_Key    Second_ Key    Log_Rec_Len   Gp_St_Pkt_No" << std::endl;
        dataFile << "             -1              0              0   " << std::setw(12) << packetCounter << std::endl;
        dataFile << "-----------------------Data---------------------" << std::endl;
    }

    dataFile.close();
}

void OdfRawFileContents::parseFileLabelData(std::bitset< 288 > dataBits,
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

void OdfRawFileContents::parseIdentifierData( std::bitset< 288 > dataBits,
                                              std::string& identifierGroupStringA,
                                              std::string& identifierGroupStringB,
                                              std::string& identifierGroupStringC )
{
    const int numBytes = 288 / 8;
    parseStringsBlockWrapper< numBytes, 8, 8, 20 >( dataBits, identifierGroupStringA, identifierGroupStringB,
                                                    identifierGroupStringC );
}

bool OdfRawFileContents::currentBlockIsHeader( std::bitset< 288 > dataBits,
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
    catch( std::runtime_error const& )
    {
        blockIsHeader =  false;
    }
    return blockIsHeader;
}

void OdfRawFileContents::readOdfFileBlock( std::istream& file, std::bitset< 36 * 8 >& dataBits )
{
    readBinaryFileBlock< 36 >( file, dataBits );
}

OdfRawFileContents::OdfRawFileContents( const std::string& odfFile ):
    fileName_( odfFile )
{
    // Open file
    std::ifstream dataFile( odfFile, std::ios_base::binary);
    if ( !dataFile.good( ) )
    {
        throw std::runtime_error( "Error when opening ODF file, file " + odfFile +  " could not be opened." );
    }

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

    // Parse file label data
    readOdfFileBlock( dataFile, currenBitBlock );
    parseFileLabelData( currenBitBlock, systemId_, programId_, spacecraftId_, fileCreationDate_,
                        fileCreationTime_, fileReferenceDate_, fileReferenceTime_ );

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
    parseIdentifierData( currenBitBlock, identifierGroupStringA_, identifierGroupStringB_, identifierGroupStringC_ );

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
    dataBlocks_.clear( );
    rampBlocks_.clear( );
    clockOffsetBlocks_.clear( );

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
                if( secondaryKey > 99 || logicalRecordLength != 1 )
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
                std::shared_ptr< OdfDataBlock > currentDataBlock = std::make_shared< OdfDataBlock >( currenBitBlock );
                if( currentDataBlock != nullptr )
                {
                    dataBlocks_.push_back( currentDataBlock );
                }
            }
            // Read ramp data
            else if ( currentBlockType == 2030 )
            {
                std::shared_ptr< OdfRampBlock > currentDataBlock = std::make_shared< OdfRampBlock >( currenBitBlock );
                if( currentDataBlock != nullptr )
                {
                    rampBlocks_[ currentRampStation ].push_back( currentDataBlock );
                }
            }
            // Clock offset data
            else if ( currentBlockType == 2040 )
            {
                std::shared_ptr< OdfClockOffsetBlock > currentDataBlock = std::make_shared< OdfClockOffsetBlock >( currenBitBlock );
                if( currentDataBlock != nullptr )
                {
                    clockOffsetBlocks_[
                            std::make_pair( currentDataBlock->getPrimaryStationId_( ),
                                            currentDataBlock->getSecondaryStationId_( ) ) ] = currentDataBlock;
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
        eofHeaderFound_ = true;
    }
    else
    {
        eofHeaderFound_ = false;
    }
}

} // namespace input_output

} // namespace tudat
