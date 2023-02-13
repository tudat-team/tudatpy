/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_READBINARYFILE_H
#define TUDAT_READBINARYFILE_H

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

//template< int NumberOfBits >
//int getSignedNBitInteger(
//        std::bitset< NumberOfBits > inputBits )
//{
//    int outputInteger = -inputBits[ NumberOfBits - 1 ] * std::pow( 2, NumberOfBits - 1 );
//
//    for( unsigned int i = 0; i < NumberOfBits - 1; i ++ )
//    {
//        outputInteger += inputBits[ i ] * std::pow( 2.0, i );
//    }
//    return outputInteger;
//}

template< unsigned int NumberOfBits >
long convertBitsetToLong(const std::bitset< NumberOfBits >& bits) {
    if ( NumberOfBits > 32 )
    {
        throw std::runtime_error( "Error when converting bitset to long: specified number of bits (" +
            std::to_string( NumberOfBits ) + "is larger than what is possible to represent with long (32).");
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

template < int NumberOfBytes >
void readBinaryFileBlock( std::istream& file,
                          std::bitset< NumberOfBytes * 8 >& dataBits )
{
    int numberOfBits = NumberOfBytes * 8;

    char dataChar [NumberOfBytes];
    file.read( (char*)dataChar, NumberOfBytes );

    if ( !file.good( ) )
    {
        throw std::runtime_error( "Error when reading data block from binary file." );
    }

    // Convert to bitset
    for ( int i = 0, bitCounter = 0; i < NumberOfBytes; ++i )
    {
        // Extract byte
        uint8_t byte = dataChar[i];
        for ( int j = 0; j < 8; ++j)
        {
            // Right shift byte to determine value of desired bit and save it to the bitset
            // Indexing of the byte and bitset starts from the right (i.e. the 0th bit is the rightmost one)
            dataBits[ numberOfBits - bitCounter - 1 ] = (byte >> ( 8 - 1 - j) ) & 1;
            ++bitCounter;
        }
    }
}

// Note: "unsignedItemFlag.at( argumentCounter )" could be replaced by "std::is_unsigned< T >::value". In that case,
// it would no longer be necessary to have unsignedItemFlag as an argument. However, that is more error-prone in case
// the signed/unsigned types aren't specified correctly.
template< unsigned int NumberBlockBits >
void parseNumericalDataBlock (std::bitset< NumberBlockBits > dataBits,
                              const std::vector< bool >& unsignedItemFlag,
                              unsigned int argumentCounter,
                              unsigned int startBitCounter )
{
    if ( startBitCounter != NumberBlockBits )
    {
        throw std::runtime_error(
                "Error when parsing binary file: block size (" + std::to_string( NumberBlockBits ) +
                " bits) and total item size (" + std::to_string(startBitCounter) + "bits ) are not consistent." );
    }
    else if ( argumentCounter != unsignedItemFlag.size() )
    {
        throw std::runtime_error(
                "Error when parsing binary file: numbers of items (" + std::to_string( argumentCounter ) +
                ") and size of unsigned flag vector (" + std::to_string( unsignedItemFlag.size() ) +
                ") are not consistent." );
    }
}

template< unsigned int NumberBlockBits, unsigned int NumberItemBits, unsigned int... NumberItemBitsN,
        typename T, typename... TN >
void parseNumericalDataBlock (std::bitset< NumberBlockBits > dataBits,
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

    parseNumericalDataBlock< NumberBlockBits, NumberItemBitsN ... >( dataBits, unsignedItemFlag, argumentCounter,
                                                                     startBitCounter, args ... );
}

template< unsigned int NumberBlockBits, unsigned int NumberItemBits, unsigned int... NumberItemBitsN,
        typename T, typename... TN >
void parseNumericalDataBlockWrapper (std::bitset< NumberBlockBits > dataBits,
                                     const std::vector< bool >& unsignedItemFlag,
                                     T& arg, TN&... args)
{
    parseNumericalDataBlock< NumberBlockBits, NumberItemBits, NumberItemBitsN ... >( dataBits, unsignedItemFlag, 0, 0,
                                                                                     arg,
                                                                                     args ... );
}

template< unsigned int NumberBlockBytes >
void parseStringsBlock (std::bitset< NumberBlockBytes * 8 > dataBits,
                        unsigned int argumentCounter,
                        unsigned int startByteCounter )
{
        if ( startByteCounter != NumberBlockBytes )
        {
            throw std::runtime_error(
                    "Error when parsing binary file: block size (" + std::to_string( NumberBlockBytes ) +
                    " bytes) and total item size (" + std::to_string( startByteCounter ) + " bytes) are not consistent." );
        }
}

template< unsigned int NumberBlockBytes, unsigned int NumberItemBytes, unsigned int... NumberItemBytesN,
        typename... TN >
void parseStringsBlock (std::bitset< NumberBlockBytes * 8 > dataBits,
                        unsigned int argumentCounter,
                        unsigned int startByteCounter,
                        std::string& arg, TN&... args )
{
    arg.resize( NumberItemBytes );

    for ( unsigned int i = 0; i < NumberItemBytes; ++i )
    {
        arg[i] = getBitsetSegment< 8, NumberBlockBytes * 8 >( dataBits, (startByteCounter + i) * 8 ).to_ulong();
    }

    ++argumentCounter;
    startByteCounter += NumberItemBytes;

    parseStringsBlock< NumberBlockBytes, NumberItemBytesN ... >( dataBits, argumentCounter,
                                                                 startByteCounter, args ...);
}

template< unsigned int NumberBlockBytes, unsigned int NumberItemBytes, unsigned int... NumberItemBytesN,
        typename... TN >
void parseStringsBlockWrapper (std::bitset< NumberBlockBytes * 8 > dataBits,
                               std::string& arg, TN&... args )
{
    parseStringsBlock< NumberBlockBytes, NumberItemBytes, NumberItemBytesN ... >(
            dataBits, 0, 0, arg, args ... );
}

} // namespace input_output

} // namespace tudat

#endif //TUDAT_READBINARYFILE_H
