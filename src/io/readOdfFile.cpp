#include "tudat/io/readOdfFile.h"

namespace tudat
{
namespace input_output
{

uint32_t convertCharactersToUnsignedInt32(
        char data[ 4 ] )
{
    uint8_t bytes[4];
    bytes[ 0 ] = data[ 0 ];
    bytes[ 1 ] = data[ 1 ];
    bytes[ 2 ] = data[ 2 ];
    bytes[ 3 ] = data[ 3 ];

    return bytes[0] << 24| (bytes[1] << 16) | (bytes[2] << 8) | (bytes[3]);
}

int32_t convertCharactersToSignedInt32(
        char data[ 4 ] )
{
    uint8_t bytes[4];
    bytes[ 0 ] = data[ 0 ];
    bytes[ 1 ] = data[ 1 ];
    bytes[ 2 ] = data[ 2 ];
    bytes[ 3 ] = data[ 3 ];

    int32_t returnData = bytes[0] << 24| (bytes[1] << 16) | (bytes[2] << 8) | (bytes[3]);
    return returnData;
}


std::shared_ptr< OdfSequentialRangeDataBlock > parseSequentialRangeData( char fileBlock[ 9 ][ 4 ], const int dopplerType )
{
    std::shared_ptr< OdfSequentialRangeDataBlock > rangeDataBlock =
            std::make_shared< OdfSequentialRangeDataBlock >( );

    std::bitset< 32 > dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 5 ] ) );

    rangeDataBlock->lowestRangingComponent = getBitsetSegment< 7, 32 >( dataBits, 0 ).to_ulong( );
    rangeDataBlock->spacecraftId = getBitsetSegment< 10, 32 >( dataBits, 7 ).to_ulong( );
    rangeDataBlock->reservedBlock = getBitsetSegment< 1, 32 >( dataBits, 17 ).to_ulong( );

    std::bitset< 14 > referenceFrequencyHighPartBitsA = getBitsetSegment< 14, 32 >( dataBits, 18 );
    dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 6 ] ) );
    std::bitset< 8 > referenceFrequencyHighPartBitsB = getBitsetSegment< 8, 32 >( dataBits, 0 );
    rangeDataBlock->referenceFrequencyHighPart = mergeBitsets< 14, 8 >(
                referenceFrequencyHighPartBitsA, referenceFrequencyHighPartBitsB ).to_ulong();
    rangeDataBlock->referenceFrequencyLowPart = getBitsetSegment< 24, 32 >( dataBits, 8 ).to_ulong( );

    dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 7 ] ) );
    rangeDataBlock->coderInPhaseTimeOffset = getSignedNBitInteger< 20 >( getBitsetSegment< 20, 32 >( dataBits, 0 ) );

    std::bitset< 12 > compositeTwoPartA = getBitsetSegment< 12, 32 >( dataBits, 20 );
    dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 8 ] ) );
    std::bitset< 10 > compositeTwoPartB = getBitsetSegment< 10, 32 >( dataBits, 0 );
    rangeDataBlock->compositeTwo = mergeBitsets< 12, 10 >(
                compositeTwoPartA, compositeTwoPartB ).to_ulong();
    rangeDataBlock->transmittingStationUplinkDelay = getBitsetSegment< 22, 32 >( dataBits, 10 ).to_ulong( );

    return rangeDataBlock;
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

    parseDataBlock< 128, 7, 10, 1, 22, 24, 20, 22, 22 >(
            dataBits, unsignedItemFlag, 0, 0,
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

    parseDataBlock< 128, 7, 10, 1, 22, 24, 20, 22, 22 >(
            dataBits, unsignedItemFlag, 0, 0,
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

    parseDataBlock< 160, 32, 10, 22, 32, 32, 3, 7, 7, 2, 6, 2, 2, 2, 1 >(
            commonDataBits, unsignedCommonItemFlag, 0, 0,
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

//            std::cout<< commonDataBlock->integerTimeTag<<" "<<
//                        commonDataBlock->fractionalTimeTag<<" "<<
//                        commonDataBlock->receivingStationDownlinkDelay<<" "<<
//                        commonDataBlock->integerObservable<<" "<<
//                        commonDataBlock->fractionalObservable<<" "<<
//                        commonDataBlock->formatId<<" "<<
//                        commonDataBlock->receivingStation<<" "<<
//                        commonDataBlock->transmittingStation<<" "<<
//                        commonDataBlock->transmittingStationNetworkId<<" "<<
//                        dataType<<" "<<
//                        commonDataBlock->downlinkBand<<" "<<
//                        commonDataBlock->uplinkBand<<" "<<
//                        commonDataBlock->referenceBand<<" "<<
//                        commonDataBlock->validity<<std::endl;

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
    parseDataBlock< 4 * 32, 32, 32, 32, 32 >(
            validDataBits, unsignedItemFlag, 0, 0,
            primaryKey,
            secondaryKey,
            logicalRecordLength,
            groupStartPacketNumber );

    // Read subsequent blocks and check if zero
    uint32_t testInteger;
    for( int i = 0; i < 5; i++ )
    {
        std::bitset< 32 > testBits = getBitsetSegment< 32, 288 >( dataBits, ( i + 4 ) * 32 );
        parseDataBlock< 32, 32 >( testBits, { true }, 0, 0, testInteger );

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
    parseDataBlock< 5 * 32, 32, 32, 32, 32, 32 >(
            numericalDataBits, unsignedItemFlag, 0, 0,
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

std::shared_ptr< OdfDopplerDataBlock > parseDopplerOrbitData( char fileBlock[ 9 ][ 4 ], const int dopplerType )
{
    std::shared_ptr< OdfDopplerDataBlock > dopplerDataBlock =
            std::make_shared< OdfDopplerDataBlock >( dopplerType );

    std::bitset< 32 > dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 5 ] ) );

    dopplerDataBlock->receiverChannel = getBitsetSegment< 7, 32 >( dataBits, 0 ).to_ulong( );
    dopplerDataBlock->spacecraftId = getBitsetSegment< 10, 32 >( dataBits, 7 ).to_ulong( );
    dopplerDataBlock->receiverExciterFlag = getBitsetSegment< 1, 32 >( dataBits, 17 ).to_ulong( );

    std::bitset< 14 > referenceFrequencyHighPartBitsA = getBitsetSegment< 14, 32 >( dataBits, 18 );
    dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 6 ] ) );
    std::bitset< 8 > referenceFrequencyHighPartBitsB = getBitsetSegment< 8, 32 >( dataBits, 0 );
    dopplerDataBlock->referenceFrequencyHighPart = mergeBitsets< 14, 8 >(
                referenceFrequencyHighPartBitsA, referenceFrequencyHighPartBitsB ).to_ulong();
    dopplerDataBlock->referenceFrequencyLowPart = getBitsetSegment< 24, 32 >( dataBits, 8 ).to_ulong( );

    dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 7 ] ) );

    dopplerDataBlock->reservedSegment = getBitsetSegment< 20, 32 >( dataBits, 0 ).to_ulong( );
    std::bitset< 12 > compressionTimeABits = getBitsetSegment< 12, 32 >( dataBits, 20 );
    dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 8 ] ) );
    std::bitset< 10 > compressionTimeBBits = getBitsetSegment< 10, 32 >( dataBits, 0 );
    dopplerDataBlock->compressionTime = mergeBitsets< 12, 10 >(
                compressionTimeABits, compressionTimeBBits ).to_ulong( );
    dopplerDataBlock->transmittingStationUplinkDelay = getBitsetSegment< 22, 32 >( dataBits, 10 ).to_ulong( );

    dopplerDataBlock->printContents();

//    std::cout<<dopplerDataBlock->referenceFrequencyHighPart<<" "<<
//               dopplerDataBlock->transmittingStationDelay<<" "<<dopplerDataBlock->receiverChannel<<" "<<
//               dopplerDataBlock->receiverExciterFlag<<" "<< dopplerDataBlock->reservedSegment <<std::endl;

    return dopplerDataBlock;
}

std::shared_ptr< OdfDataBlock > parseOrbitData( char fileBlock[ 9 ][ 4 ] )
{
    std::shared_ptr< OdfCommonDataBlock > commonDataBlock =
            std::make_shared< OdfCommonDataBlock >( );

    commonDataBlock->integerTimeTag = convertCharactersToUnsignedInt32( fileBlock[ 0 ] );

    std::bitset< 32 > dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 1 ] ) );
    commonDataBlock->fractionalTimeTag = getBitsetSegment< 10, 32 >(
                dataBits, 0 ).to_ulong( );
    commonDataBlock->receivingStationDownlinkDelay = getBitsetSegment< 22, 32 >(
                dataBits, 10 ).to_ulong( );

    commonDataBlock->integerObservable = convertCharactersToSignedInt32( fileBlock[ 2 ] );
    commonDataBlock->fractionalObservable = convertCharactersToSignedInt32( fileBlock[ 3 ] );


    dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 4 ] ) );


    commonDataBlock->formatId = getBitsetSegment< 3, 32 >( dataBits, 0 ).to_ulong( );
    commonDataBlock->receivingStation = getBitsetSegment< 7, 32 >( dataBits, 3 ).to_ulong( );
    commonDataBlock->transmittingStation = getBitsetSegment< 7, 32 >( dataBits, 10 ).to_ulong( );

    commonDataBlock->transmittingStationNetworkId = getBitsetSegment< 2, 32 >( dataBits, 17 ).to_ulong( );
    int dataType =  getBitsetSegment< 6, 32 >( dataBits, 19 ).to_ulong( );

    commonDataBlock->downlinkBand = getBitsetSegment< 2, 32 >( dataBits, 25 ).to_ulong( );
    commonDataBlock->uplinkBand = getBitsetSegment< 2, 32 >( dataBits, 27 ).to_ulong( );
    commonDataBlock->referenceBand = getBitsetSegment< 2, 32 >( dataBits, 29 ).to_ulong( );
    commonDataBlock->validity = getBitsetSegment< 1, 32 >( dataBits, 31 ).to_ulong( );

    std::shared_ptr< OdfDataBlock > dataBlock = std::make_shared< OdfDataBlock >( );


    dataBlock->commonDataBlock = commonDataBlock;
    if(  dataType == 11 || dataType == 12 || dataType == 13 )
    {
        if( dataType != 11 )
        {
//            std::cout<<commonDataBlock->integerObservable<<" "<<
//                       commonDataBlock->downlinkBand<<" "<<
//                       commonDataBlock->uplinkBand<<" "<<
//                       commonDataBlock->referenceBand<<" "<<
//                       commonDataBlock->validity<<" "<<
//                       commonDataBlock->formatId<<" "<<
//                       commonDataBlock->transmittingStationNetworkId<<" "<<
//                       commonDataBlock->receivingStation<<" "<<
//                       commonDataBlock->transmittingStation<<" "<<
//                       commonDataBlock->receivingStationDownlinkDelay<<" ";
        }
        std::bitset< 32 > dataBits1 = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 5 ] ) );
        std::bitset< 32 > dataBits2 = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 6 ] ) );
        std::bitset< 32 > dataBits3 = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 7 ] ) );
        std::bitset< 32 > dataBits4 = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 8 ] ) );
        std::bitset< 128 > dataBits5 = mergeBitsets<96, 32>( mergeBitsets<64, 32>(
                mergeBitsets<32, 32>(dataBits1, dataBits2), dataBits3 ), dataBits4 );

        dataBlock->observableSpecificDataBlock = parseDopplerOrbitData( dataBits5, dataType );

//        dataBlock->observableSpecificDataBlock = parseDopplerOrbitData( fileBlock, dataType );
    }
//    else if(  dataType == 37 )
//    {
//        dataBlock->observableSpecificDataBlock = parseSequentialRangeData( fileBlock, dataType );
//    }
    else
    {
        //std::cerr<<"Data type "<<dataType<<std::endl;
        dataBlock = nullptr;
        //throw std::runtime_error( "Error, ODF data type " + std::to_string( dataType ) + " not recognized" );

    }
    return dataBlock;
}

OdfRampBlock parseRampData( char fileBlock[ 9 ][ 4 ] )
{

    OdfRampBlock rampBlock;
    rampBlock.integerRampStartTime = convertCharactersToUnsignedInt32( fileBlock[ 0 ] );
    rampBlock.fractionalRampStartTime = convertCharactersToUnsignedInt32( fileBlock[ 1 ] );

    rampBlock.integerRampRate = convertCharactersToSignedInt32( fileBlock[ 2 ] );
    rampBlock.fractionalRampRate = convertCharactersToSignedInt32( fileBlock[ 3 ] );

    std::bitset< 32 > dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 4 ] ) );
    rampBlock.transmittingStationId = getBitsetSegment< 10, 32 >( dataBits, 22 ).to_ulong( );
    rampBlock.integerRampStartFrequency = getBitsetSegment< 22, 32 >( dataBits, 0 ).to_ulong( );

    rampBlock.integerRampStartFrequencyModulo = convertCharactersToSignedInt32( fileBlock[ 5 ] );
    rampBlock.fractionalRampStartFrequency = convertCharactersToSignedInt32( fileBlock[ 6 ] );

    rampBlock.integerRampEndTime = convertCharactersToSignedInt32( fileBlock[ 7 ] );
    rampBlock.fractionalRampEndTime = convertCharactersToSignedInt32( fileBlock[ 8 ] );


    return rampBlock;
}

void parseFileLabel( char fileBlock[ 9 ][ 4 ],
std::string& systemId, std::string& programId, std::string& fileCreationDate, std::string& fileCreationTme,
uint32_t& spacecraftIdNumber, uint32_t& fileReferenceDate, uint32_t& fileReferenceTime )
{
    systemId.resize( 8 );
    programId.resize( 8 );

    fileCreationDate.resize( 4 );
    fileCreationTme.resize( 4 );

    // Read strong contents of data block
    for( unsigned int i = 0; i < 4; i++ )
    {
        systemId[ i ] = fileBlock[ 0 ][ i ];
        systemId[ i + 4 ] = fileBlock[ 1 ][ i ];

        programId[ i ] = fileBlock[ 2 ][ i ];
        programId[ i + 4 ] = fileBlock[ 3 ][ i ];

        fileCreationDate[ i ] = fileBlock[ 5 ][ i ];
        fileCreationTme[ i ] = fileBlock[ 6 ][ i ];
    }

    // Retrieve numeral contents of data block
    spacecraftIdNumber = convertCharactersToUnsignedInt32( fileBlock[ 4 ] );
    fileReferenceDate = convertCharactersToUnsignedInt32( fileBlock[ 7 ] );
    fileReferenceTime = convertCharactersToUnsignedInt32( fileBlock[ 8 ] );
}

void parseHeader( char fileBlock[ 9 ][ 4 ],
int32_t& primaryKey,
uint32_t& secondaryKey,
uint32_t& logicalrecordLength,
uint32_t& groupStartPacketNumber )
{
    // Read data from header file
    primaryKey = convertCharactersToSignedInt32( fileBlock[ 0 ] );
    secondaryKey = convertCharactersToUnsignedInt32( fileBlock[ 1 ] );
    logicalrecordLength = convertCharactersToUnsignedInt32( fileBlock[ 2 ] );
    groupStartPacketNumber = convertCharactersToUnsignedInt32( fileBlock[ 3 ] );

    // Check if subsequent data contains only zeroes
    uint32_t testInteger;
    for( int i = 0; i < 5; i++ )
    {
        testInteger = convertCharactersToUnsignedInt32( fileBlock[ i + 4 ] );
        if( testInteger != 0 )
        {
            throw std::runtime_error( "Error when reading ODF file, header file inconsistent" );
        }
    }
}


int currentBlockIsHeader(  char fileBlock[ 9 ][ 4 ], unsigned int& secondaryKeyInt )
{
    int headerType = -2;

    int32_t primaryKey;
    uint32_t secondaryKey;
    uint32_t logicalrecordLength;
    uint32_t groupStartPacketNumber;

    // Try reading current block as a header, if an exception is caught, block is not a header.
    try
    {
        parseHeader( fileBlock, primaryKey, secondaryKey, logicalrecordLength, groupStartPacketNumber );

        // If block is a header, retrieve secondary key and header type
        secondaryKeyInt = secondaryKey;
        if( primaryKey ==  101 )
        {
            headerType = 1;
        }
        else if( primaryKey == 107 )
        {
            headerType = 2;
        }
        else if( primaryKey == 109 )
        {
            headerType = 3;
        }
        else if( primaryKey == 2030 )
        {
            headerType = 4;
        }
        else if( primaryKey == 2040 )
        {
            headerType = 5;
        }
        else if( primaryKey == -1 )
        {
            headerType = -1;
        }
    }
    catch( std::runtime_error )
    {
        headerType = 0;
    }
    return headerType;
}


void parseIdentifierGroup( char fileBlock[ 9 ][ 4 ],
    std::string& identifierGroupStringA, std::string&identifierGroupStringB, std::string& identifierGroupStringC )
{
    // Read three strings from data blocks
    identifierGroupStringA.resize( 8 );
    identifierGroupStringB.resize( 8 );
    identifierGroupStringC.resize( 20 );

    for( unsigned int i = 0; i < 4; i++ )
    {
        identifierGroupStringA[ i ] = fileBlock[ 0 ][ i ];
        identifierGroupStringA[ i + 4 ] = fileBlock[ 1 ][ i ];

        identifierGroupStringB[ i ] = fileBlock[ 2 ][ i ];
        identifierGroupStringB[ i + 4 ] = fileBlock[ 3 ][ i ];

        for( unsigned int j = 0; j < 5; j++ )
        {
            identifierGroupStringC[ i + j * 4 ] = fileBlock[ 4 + j ][ i ];
        }
    }
}

//! Function to read a single 36 byte block from ODF file
void readOdfFileBlock(
        char fileBlock[ 9 ][ 4 ],
std::istream& file )
{
    for( unsigned int i = 0; i < 9; i++ )
    {
        file.read( (char*)fileBlock[ i ], 4 );
    }

    std::cout << fileBlock << std::endl;
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

    // Read first line
//    char currentFileBlock[ 9 ][ 4 ];
//    readOdfFileBlock( currentFileBlock, dataFile );

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

//    std::cout<< odfFileContents->systemId << " " << odfFileContents->programId << " " << odfFileContents->spacecraftId <<
//        " " << odfFileContents->fileCreationDate << " " << odfFileContents->fileCreationTime << " " <<
//        odfFileContents->fileReferenceDate << " " << odfFileContents->fileReferenceTime << std::endl;

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

    // Define output of data blocks and ramp blocks
    std::vector< std::shared_ptr< OdfDataBlock > > unsortedOdfDataBlocks;
    std::map< int, std::vector< OdfRampBlock > > odfRampBlocks;

    bool continueFileRead = true;
    int currentRampStation = -1;
    int currentBlockType = 109;

    // Read file until end
    while( continueFileRead && !dataFile.eof( ) )
    {
        // Read current block
        readOdfFileBlock( dataFile, currentFileBlock );

        // Check if block is a header file
        bool blockIsHeader = currentBlockIsHeader(
                currentFileBlock, primaryKey, secondaryKey, logicalRecordLength, groupStartPacketNumber );

        // If block is header
        if ( blockIsHeader )
        {
            throw std::runtime_error( "DUMMY Error when reading ODF group, invalid block type." );
            // Ramp group header
            if ( primaryKey == 2030 )
            {
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
                    unsortedOdfDataBlocks.push_back( currentDataBlock );
                }
            }
            // Read ramp data
            else if ( currentBlockType == 2030 )
            {

            }
            else
            {
                throw std::runtime_error( "Error when reading ODF group, invalid block type." );
            }
        }

//        if( blockIsHeader == 0 && currentBlockType == 3 )
//        {
//            std::shared_ptr< OdfDataBlock > currentDataBlock = parseOrbitData( currentFileBlock );
//            if( currentDataBlock != NULL )
//            {
//                unsortedOdfDataBlocks.push_back( currentDataBlock );
//            }
//        }
//        else if( blockIsHeader == 0 && currentBlockType == 4 )
//        {
//            odfRampBlocks[ currentRampStation ].push_back( parseRampData( currentFileBlock ) );
//        }
//        // If block is a header, set current data block type accordingly.
//        else if( blockIsHeader != 0 )
//        {
//            currentBlockType = blockIsHeader;
//            if( currentBlockType == 4 )
//            {
//                currentRampStation = secondaryKey;
//            }
//            if( currentBlockType == -1 )
//            {
//                continueFileRead = 0;
//            }
//        }
//        else
//        {
//            throw std::runtime_error( "Error, did not recognized header ODF block" );
//        }
    }

    // Set file contents
    odfFileContents->dataBlocks = unsortedOdfDataBlocks;
    odfFileContents->odfRampBlocks = odfRampBlocks;


//    if( continueFileRead == 1 )
//    {
//        odfFileContents->eofHeaderFound = false;
//    }
//    else
//    {
//        odfFileContents->eofHeaderFound = true;
//    }

    return odfFileContents;
}

} // namespace input_output

} // namespace tudat
