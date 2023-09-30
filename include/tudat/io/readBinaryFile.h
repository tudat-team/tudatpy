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

//! Function to concatenate two bitsets.
/*!
 * @tparam FirstInputSize Size of the first bitset.
 * @tparam SecondInputSize Size of the second bitset.
 * @param firstInput Most significant bitset.
 * @param secondInput Least significant bitset.
 * @return Concatenated bitset.
 */
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

//! Function to extract a segment from a bitset.
/*!
 * @tparam OutputBits Desired size of outputted bitset.
 * @tparam InputBits Size of inputted bitset.
 * @param inputBits Bitset from which a segment will be extracted.
 * @param startIndex Index from which to start extracting the desired bitset. 0 corresponds to the leftmost bit (i.e.
 * most significant bit)
 * @return Segment of inputBits.
 */
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

//! Function to convert an arbitrarily sized bitset to a long (i.e. binary to signed decimal number).
/*!
 * @tparam NumberOfBits Size of the provided bitset
 * @param bits Bitset
 * @return Signed integer representation of the bitset.
 */
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

//! Function to a specified number of bytes from a file into a bitset.
/*!
 * Function to a specified number of bytes from a file into a bitset. The read bitset is returned by reference.
 * @tparam NumberOfBytes Number of bytes to read from the file.
 * @param file Input files.
 * @param dataBits Bitset read from the file.
 */
template < int NumberOfBytes >
void readBinaryFileBlock( std::istream& file,
                          std::bitset< NumberOfBytes * 8 >& dataBits )
{
    int numberOfBits = NumberOfBytes * 8;

    char dataChar [NumberOfBytes];
    file.read( static_cast< char* >( dataChar ), NumberOfBytes );

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
            dataBits[ numberOfBits - bitCounter - 1 ] = ( byte >> ( 8 - 1 - j ) ) & 1;
            ++bitCounter;
        }
    }
}

//! Function to extract an arbitrary number of signed and unsigned integers from a bitset.
/*!
 * Function to extract an arbitrary number of integers from a bitset. Overload for 0 arguments, simply executes some
 * error checking. Should not be called directly, but instead via parseNumericalDataBlockWrapper.
 *
 * @tparam NumberBlockBits Size of the inputted bitset.
 * @param dataBits Bitset from which to extract data
 * @param unsignedItemFlag Vector indicating whether each number should be extracted as signed or unsigned
 * @param argumentCounter Counter for the number of integers already extracted
 * @param startBitCounter Start bit of the current integer
 * provided NumberItemBits, NumberItemBitsN.
 */
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

//! Function to extract an arbitrary number of signed and unsigned integers from a bitset.
/*!
 * Function to extract an arbitrary number of signed and unsigned integers from a bitset. The extracted numbers are
 * returned by reference. Overload for >= 1 arguments. Should not be called directly, but instead via
 * parseNumericalDataBlockWrapper.
 *
 * Note: "unsignedItemFlag.at( argumentCounter )" could be replaced by "std::is_unsigned< T >::value". In that case,
 * it would no longer be necessary to have unsignedItemFlag as an argument. However, that is more error-prone in case
 * the signed/unsigned types aren't specified correctly.
 *
 * @tparam NumberBlockBits Size of the inputted bitset.
 * @tparam NumberItemBits, NumberItemBitsN Number of bits to convert into each integer.
 * @tparam T, TN Types of the extracted integers. The type doesn't influence the conversion from binary to decimal.
 * @param dataBits Bitset from which to extract data
 * @param unsignedItemFlag Vector indicating whether each number should be extracted as signed or unsigned
 * @param argumentCounter Counter for the number of integers already extracted
 * @param startBitCounter Start bit of the current integer
 * @param arg, args Extracted signed/unsigned integers. Their number should coincide with the number of
 * provided NumberItemBits, NumberItemBitsN.
 */
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

//! Function to extract an arbitrary number of signed and unsigned integers from a bitset.
/*!
 * Function to extract an arbitrary number of integers from a bitset. The extracted numbers are returned by reference.
 *
 * @tparam NumberBlockBits Size of the inputted bitset.
 * @tparam NumberItemBits, NumberItemBitsN Number of bits to convert into each integer.
 * @tparam T, TN Types of the extracted integers. The type doesn't influence the conversion from binary to decimal.
 * @param dataBits Bitset from which to extract data
 * @param unsignedItemFlag Vector indicating whether each number should be extracted as signed or unsigned
 * @param arg, args Extracted signed/unsigned integers. Their number should coincide with the number of
 * provided NumberItemBits, NumberItemBitsN.
 */
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

//! Function to extract an arbitrary number of strings from a bitset.
/*!
 * Function to extract an arbitrary number of strings from a bitset. Overload for 0 arguments, simply executes some
 * error checking. Should not be called directly, but instead via parseStringsBlockWrapper.
 *
 * @tparam NumberBlockBytes Size of the inputted bitset in bytes.
 * @param dataBits Bitset from which to extract data
 * @param argumentCounter Counter for the number of strings already extracted
 * @param startByteCounter Start byte of the current string
 */
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

//! Function to extract an arbitrary number of strings from a bitset.
/*!
 * Function to extract an arbitrary number of strings from a bitset. The extracted strings are returned by reference.
 * Overload for >=1 arguments. Should not be called directly, but instead via parseStringsBlockWrapper.
 *
 * @tparam NumberBlockBytes Size of the inputted bitset in bytes.
 * @tparam NumberItemBytes, NumberItemBytesN Number of bytes to convert into each string.
 * @tparam TN Type of args. Has to be std::string (enforced via the variadic template "recursion")
 * @param dataBits Bitset from which to extract data
 * @param argumentCounter Counter for the number of strings already extracted
 * @param startByteCounter Start byte of the current string
 * @param arg, args Extracted strings. Their number should coincide with the number of the items
 * provided in NumberItemBytes, NumberItemBytesN.
 */
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

//! Function to extract an arbitrary number of strings from a bitset.
/*!
 * Function to extract an arbitrary number of strings from a bitset. The extracted strings are returned by reference.
 *
 * @tparam NumberBlockBytes Size of the inputted bitset in bytes.
 * @tparam NumberItemBytes, NumberItemBytesN Number of bytes to convert into each string.
 * @tparam TN Type of args. Has to be std::string (enforced via the variadic template "recursion")
 * @param dataBits Bitset from which to extract data
 * @param arg, args Extracted strings. Their number should coincide with the number of the items
 * provided in NumberItemBytes, NumberItemBytesN.
 */
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
