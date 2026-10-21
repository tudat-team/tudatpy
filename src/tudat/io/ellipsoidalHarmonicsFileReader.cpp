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
#include <stdexcept>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/math/special_functions/binomial.hpp>

namespace tudat
{

namespace input_output
{

using namespace tudat::basic_mathematics;

namespace
{

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

//! Allocates one square matrix per group per degree, sized by the number of terms in the group.
/*!
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

//! Reads one of the two files that hold a line per term, "degree orderIndex label values...".
/*!
 * \param fileName File to read.
 * \param maximumDegree Highest degree to read.
 * \param storeAsColumn If true each term's values become a column of its group's matrix, which
 *      is how the Lame polynomials are laid out. If false they become a row, as the transform
 *      factors are.
 * \param convertPolynomial If true the values are converted from powers of (1 - x) to powers
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

    std::ifstream stream( fileName );
    if( !stream.good( ) )
    {
        throw std::runtime_error( "Error when reading ellipsoidal harmonics file, could not open " + fileName );
    }
    std::string line;
    int lineNumber = 0;

    while( std::getline( stream, line ) )
    {
        lineNumber++;
        std::vector< std::string > fields;
        boost::algorithm::split(
                fields, boost::algorithm::trim_copy( line ), boost::algorithm::is_space( ), boost::algorithm::token_compress_on );
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

        Eigen::VectorXd values( fields.size( ) - 3 );
        for( int i = 0; i < values.size( ); i++ )
        {
            values( i ) = parseFortranDouble( fields.at( i + 3 ) );
        }

        // The group label, the number of values and the term itself must match what this reader expects.
        Eigen::MatrixXd& matrix = matrices.getGroup( degree, group );
        if( fields.at( 2 ).at( 0 ) != getGroupLabel( group ) || values.size( ) != matrix.rows( ) ||
            isTermRead[ degree ][ fileOrderIndex - 1 ] )
        {
            throw std::runtime_error( "Error when reading ellipsoidal harmonics file " + fileName + ", line " +
                                      std::to_string( lineNumber ) + " should hold the only entry for degree " + std::to_string( degree ) +
                                      " order index " + std::to_string( fileOrderIndex ) + ", with group label " +
                                      std::string( 1, getGroupLabel( group ) ) + " and " + std::to_string( matrix.rows( ) ) + " values." );
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
    const std::string converted = boost::algorithm::replace_all_copy( boost::algorithm::to_upper_copy( text ), "D", "E" );

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
            sum += coefficients( j ) * boost::math::binomial_coefficient< double >( j, i );
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

    std::ifstream stream( fileName );
    if( !stream.good( ) )
    {
        throw std::runtime_error( "Error when reading ellipsoidal harmonics file, could not open " + fileName );
    }
    std::string line;

    while( std::getline( stream, line ) )
    {
        std::vector< std::string > fields;
        boost::algorithm::split(
                fields, boost::algorithm::trim_copy( line ), boost::algorithm::is_space( ), boost::algorithm::token_compress_on );
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

        if( isOrderRead[ degree ][ order ] )
        {
            throw std::runtime_error( "Error when reading ellipsoidal harmonic coefficients from " + fileName + ", degree " +
                                      std::to_string( degree ) + " order " + std::to_string( order ) + " appears twice." );
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
