/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *          Garmier, R., Barriot, J.-P. (2001), Celestial Mechanics and Dynamical Astronomy, 79(4), 235-275.
 *          Hu, X., Jekeli, C. (2015), Journal of Geodesy, 89(2), 159-177.
 *          See tudat/math/basic/ellipsoidalHarmonics.h for the full citations.
 */

#include "tudat/io/ellipsoidalHarmonicsFileReader.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace tudat
{

namespace input_output
{

using namespace tudat::basic_mathematics;

namespace
{

//! Opens a file, or throws naming it.
std::ifstream openOrThrow( const std::string& fileName )
{
    std::ifstream stream( fileName );
    if( !stream.good( ) )
    {
        throw std::runtime_error( "Error when reading ellipsoidal harmonics file, could not open " + fileName );
    }
    return stream;
}

//! Splits a line on whitespace.
std::vector< std::string > splitLine( const std::string& line )
{
    std::vector< std::string > fields;
    std::istringstream stream( line );
    std::string field;
    while( stream >> field )
    {
        fields.push_back( field );
    }
    return fields;
}

//! Binomial coefficient, for the small values the polynomial conversion needs.
double binomialCoefficient( const int upper, const int lower )
{
    if( lower < 0 || lower > upper )
    {
        return 0.0;
    }

    double value = 1.0;
    for( int i = 0; i < lower; i++ )
    {
        value = value * ( upper - i ) / ( i + 1 );
    }
    return value;
}

//! Letter used in the Lame and transform files for each term group.
char getGroupLabel( const EllipsoidalHarmonicTermGroup group )
{
    switch( group )
    {
        case ellipsoidal_harmonic_k_terms:
            return 'K';
        case ellipsoidal_harmonic_l_terms:
            return 'L';
        case ellipsoidal_harmonic_m_terms:
            return 'M';
        default:
            return 'N';
    }
}

/*! Allocates one square matrix per group per degree, sized by the number of terms in the group.
 *
 * The Lame polynomials and the transform factors share this shape: a group with t terms has a
 * t by t matrix, holding t polynomials of t coefficients, or a t by t mixing matrix.
 */
EllipsoidalHarmonicMatrices createSquareMatrixSet( const int maximumDegree )
{
    EllipsoidalHarmonicMatrices matrices;
    matrices.kTerms.resize( maximumDegree + 1 );
    matrices.lTerms.resize( maximumDegree + 1 );
    matrices.mTerms.resize( maximumDegree + 1 );
    matrices.nTerms.resize( maximumDegree + 1 );

    for( int degree = 0; degree <= maximumDegree; degree++ )
    {
        for( int groupIndex = 0; groupIndex < 4; groupIndex++ )
        {
            const EllipsoidalHarmonicTermGroup group = static_cast< EllipsoidalHarmonicTermGroup >( groupIndex );
            const int numberOfTerms = getNumberOfTermsInGroup( degree, group );
            matrices.getGroup( degree, group ) = Eigen::MatrixXd::Zero( numberOfTerms, numberOfTerms );
        }
    }

    return matrices;
}

/*! Reads one of the two files that hold a line per term, "degree orderIndex label values...".
 *
 * @param fileName File to read.
 * @param maximumDegree Highest degree to read.
 * @param storeAsColumn If true each term's values become a column of its group's matrix, which
 *      is how the Lame polynomials are laid out; if false they become a row, as the transform
 *      factors are.
 * @param convertPolynomial If true the values are converted from powers of (1 - x) to powers
 *      of x, which the Lame polynomials need and the transform factors do not.
 */
EllipsoidalHarmonicMatrices readPerTermFile( const std::string& fileName,
                                             const int maximumDegree,
                                             const bool storeAsColumn,
                                             const bool convertPolynomial )
{
    EllipsoidalHarmonicMatrices matrices = createSquareMatrixSet( maximumDegree );

    // Which terms the file has supplied, indexed by degree and file order index. A term left
    // out would stay zero, and a zero Lame polynomial turns the whole field into NaN.
    std::vector< std::vector< bool > > isTermRead( maximumDegree + 1 );
    for( int degree = 0; degree <= maximumDegree; degree++ )
    {
        isTermRead[ degree ].assign( 2 * degree + 1, false );
    }

    std::ifstream stream = openOrThrow( fileName );
    std::string line;
    int lineNumber = 0;

    while( std::getline( stream, line ) )
    {
        lineNumber++;
        const std::vector< std::string > fields = splitLine( line );
        if( fields.size( ) < 4 )
        {
            continue;
        }

        const int degree = std::stoi( fields.at( 0 ) );
        const int fileOrderIndex = std::stoi( fields.at( 1 ) );

        if( degree > maximumDegree || degree < 0 || fileOrderIndex < 1 || fileOrderIndex > 2 * degree + 1 )
        {
            continue;
        }

        EllipsoidalHarmonicTermGroup group;
        int indexInGroup;
        getTermGroupAndIndex( degree, fileOrderIndex, group, indexInGroup );

        // The file also labels the group; disagreement means the file layout is not what this
        // reader assumes, which would otherwise show up only as a wrong field.
        if( !fields.at( 2 ).empty( ) && fields.at( 2 ).at( 0 ) != getGroupLabel( group ) )
        {
            throw std::runtime_error( "Error when reading ellipsoidal harmonics file " + fileName + ", line " +
                                      std::to_string( lineNumber ) + ": degree " + std::to_string( degree ) + " order index " +
                                      std::to_string( fileOrderIndex ) + " should be group " + std::string( 1, getGroupLabel( group ) ) +
                                      " but the file says " + fields.at( 2 ) + "." );
        }

        Eigen::VectorXd values( fields.size( ) - 3 );
        for( int i = 0; i < values.size( ); i++ )
        {
            values( i ) = parseFortranDouble( fields.at( i + 3 ) );
        }

        Eigen::MatrixXd& matrix = matrices.getGroup( degree, group );
        if( values.size( ) != matrix.rows( ) )
        {
            throw std::runtime_error( "Error when reading ellipsoidal harmonics file " + fileName + ", line " +
                                      std::to_string( lineNumber ) + ": expected " + std::to_string( matrix.rows( ) ) +
                                      " values for degree " + std::to_string( degree ) + " group " +
                                      std::string( 1, getGroupLabel( group ) ) + " but found " + std::to_string( values.size( ) ) + "." );
        }

        if( convertPolynomial )
        {
            values = convertPolynomialToPowersOfArgument( values );
        }

        if( storeAsColumn )
        {
            matrix.col( indexInGroup ) = values;
        }
        else
        {
            matrix.row( indexInGroup ) = values.transpose( );
        }
        isTermRead[ degree ][ fileOrderIndex - 1 ] = true;
    }

    for( int degree = 0; degree <= maximumDegree; degree++ )
    {
        for( int fileOrderIndex = 1; fileOrderIndex <= 2 * degree + 1; fileOrderIndex++ )
        {
            if( !isTermRead[ degree ][ fileOrderIndex - 1 ] )
            {
                throw std::runtime_error( "Error when reading ellipsoidal harmonics file " + fileName + ", no entry for degree " +
                                          std::to_string( degree ) + " order index " + std::to_string( fileOrderIndex ) +
                                          ". The file may not reach the requested maximum degree of " + std::to_string( maximumDegree ) +
                                          "." );
            }
        }
    }

    return matrices;
}

}  // namespace

double parseFortranDouble( const std::string& text )
{
    std::string converted = text;
    for( char& character : converted )
    {
        if( character == 'D' || character == 'd' )
        {
            character = 'E';
        }
    }

    try
    {
        return std::stod( converted );
    }
    catch( const std::exception& )
    {
        throw std::runtime_error( "Error when reading ellipsoidal harmonics file, could not read '" + text + "' as a number." );
    }
}

Eigen::VectorXd convertPolynomialToPowersOfArgument( const Eigen::VectorXd& coefficients )
{
    // Expanding sum_j b_j (1-x)^j and collecting powers of x gives
    // c_i = (-1)^i sum_{j >= i} b_j * binomial( j, i ).
    const int numberOfCoefficients = coefficients.size( );
    Eigen::VectorXd converted = Eigen::VectorXd::Zero( numberOfCoefficients );

    for( int i = 0; i < numberOfCoefficients; i++ )
    {
        double sum = 0.0;
        for( int j = i; j < numberOfCoefficients; j++ )
        {
            sum += coefficients( j ) * binomialCoefficient( j, i );
        }
        converted( i ) = ( i % 2 == 0 ) ? sum : -sum;
    }

    return converted;
}

EllipsoidalHarmonicCoefficients readEllipsoidalHarmonicCoefficients( const std::string& fileName, const int maximumDegree )
{
    // Values are first collected per degree in the order the file implies: cosine terms of
    // descending order, the zero-order term, then sine terms of ascending order. The four
    // groups then take consecutive slices of that.
    std::vector< Eigen::VectorXd > perDegree( maximumDegree + 1 );
    std::vector< std::vector< bool > > isOrderRead( maximumDegree + 1 );
    for( int degree = 0; degree <= maximumDegree; degree++ )
    {
        perDegree[ degree ] = Eigen::VectorXd::Zero( 2 * degree + 1 );
        isOrderRead[ degree ].assign( degree + 1, false );
    }

    std::ifstream stream = openOrThrow( fileName );
    std::string line;

    while( std::getline( stream, line ) )
    {
        const std::vector< std::string > fields = splitLine( line );
        if( fields.size( ) != 4 || fields.at( 0 ).at( 0 ) == '#' )
        {
            continue;
        }

        const int degree = std::stoi( fields.at( 0 ) );
        const int order = std::stoi( fields.at( 1 ) );

        if( degree < 0 || degree > maximumDegree || order < 0 || order > degree )
        {
            continue;
        }

        const double cosineTerm = parseFortranDouble( fields.at( 2 ) );
        const double sineTerm = parseFortranDouble( fields.at( 3 ) );

        perDegree[ degree ]( degree - order ) = cosineTerm;
        if( order > 0 )
        {
            perDegree[ degree ]( degree + order ) = sineTerm;
        }
        isOrderRead[ degree ][ order ] = true;
    }

    // A missing line would silently become a zero coefficient.
    for( int degree = 0; degree <= maximumDegree; degree++ )
    {
        for( int order = 0; order <= degree; order++ )
        {
            if( !isOrderRead[ degree ][ order ] )
            {
                throw std::runtime_error( "Error when reading ellipsoidal harmonic coefficients from " + fileName +
                                          ", no entry for degree " + std::to_string( degree ) + " order " + std::to_string( order ) +
                                          ". The file may not reach the requested maximum degree of " + std::to_string( maximumDegree ) +
                                          "." );
            }
        }
    }

    EllipsoidalHarmonicCoefficients coefficients = createEllipsoidalHarmonicCoefficients( maximumDegree );
    for( int degree = 0; degree <= maximumDegree; degree++ )
    {
        int offset = 0;
        for( int groupIndex = 0; groupIndex < 4; groupIndex++ )
        {
            const EllipsoidalHarmonicTermGroup group = static_cast< EllipsoidalHarmonicTermGroup >( groupIndex );
            const int numberOfTerms = getNumberOfTermsInGroup( degree, group );
            coefficients.getGroup( degree, group ) = perDegree[ degree ].segment( offset, numberOfTerms );
            offset += numberOfTerms;
        }
    }

    return coefficients;
}

EllipsoidalHarmonicMatrices readLameFunctionCoefficients( const std::string& fileName, const int maximumDegree )
{
    return readPerTermFile( fileName, maximumDegree, true, true );
}

EllipsoidalHarmonicMatrices readEllipsoidalHarmonicTransformFactors( const std::string& fileName, const int maximumDegree )
{
    return readPerTermFile( fileName, maximumDegree, false, false );
}

}  // namespace input_output

}  // namespace tudat
