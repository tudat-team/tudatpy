/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120926    E. Dekens         File created.
 *      121218    S. Billemont      Added output fuctions to display Legendre polynomial data,
 *                                  for debugging.
 *
 *    References
 *
 *    Notes
 *
 */

#include <sstream>
#include <stdexcept>

#include <boost/exception/all.hpp>

#include "Tudat/Mathematics/BasicMathematics/legendrePolynomials.h"

namespace tudat
{
namespace basic_mathematics
{

// Define maximum size of Legendre polynomials back-end cache.
#ifndef MAXIMUM_CACHE_ENTRIES
#define MAXIMUM_CACHE_ENTRIES 12000
#endif

//! Compute unnormalized associated Legendre polynomial.
double computeLegendrePolynomial( const int degree,
                                  const int order,
                                  const double polynomialParameter )
{
    // If degree or order is negative...
    if ( degree < 0 || order < 0 )
    {
        // Set error message.
        std::stringstream errorMessage;
        errorMessage << "Error: the Legendre polynomial of = " << degree << " and order = "
                     << order << " is undefined." << std::endl;

        // Throw a run-time error.
        boost::throw_exception( boost::enable_error_info( std::runtime_error(
               errorMessage.str( ) ) ) );
    }

    // Else if order is greater than degree...
    else if ( order > degree && degree >= 0 )
    {
        // Return zero.
        return 0.0;
    }

    // Else if order and degree are lower than 2...
    else if ( degree <= 1 && order <= 1 )
    {
        // Compute polynomial explicitly.
        return computeLegendrePolynomialExplicit( degree, order, polynomialParameter );
    }

    // Else if degree and order are sectoral...
    else if ( degree == order )
    {
        // Obtain polynomial of degree one and order one.
        const double degreeOneOrderOnePolynomial = legendreCache.getOrElseUpdate(
                    1, 1, polynomialParameter, &computeLegendrePolynomial );

        // Obtain prior sectoral polynomial.
        const double priorSectoralPolynomial = legendreCache.getOrElseUpdate(
                    degree - 1, order - 1, polynomialParameter, &computeLegendrePolynomial );

        // Compute polynomial.
        return computeLegendrePolynomialDiagonal(
                    degree, degreeOneOrderOnePolynomial, priorSectoralPolynomial );
    }

    // Else degree and order are zonal/tessoral...
    else
    {
        // Obtain prior degree polynomial.
        const double oneDegreePriorPolynomial = legendreCache.getOrElseUpdate(
                    degree - 1, order, polynomialParameter, &computeLegendrePolynomial );

        // Obtain two degrees prior polynomial.
        const double twoDegreesPriorPolynomial = legendreCache.getOrElseUpdate(
                    degree - 2, order, polynomialParameter, &computeLegendrePolynomial );

        // Compute polynomial.
        return computeLegendrePolynomialVertical( degree,
                                                  order,
                                                  polynomialParameter,
                                                  oneDegreePriorPolynomial,
                                                  twoDegreesPriorPolynomial );
    }
}

//! Compute geodesy-normalized associated Legendre polynomial.
double computeGeodesyLegendrePolynomial( const int degree,
                                         const int order,
                                         const double polynomialParameter )
{
    // If degree or order is negative...
    if ( degree < 0 || order < 0 )
    {
        // Set error message.
        std::stringstream errorMessage;
        errorMessage << "Error: the Legendre polynomial of = " << degree << " and order = "
                     << order << " is undefined." << std::endl;

        // Throw a run-time error.
        boost::throw_exception( boost::enable_error_info( std::runtime_error(
               errorMessage.str( ) ) ) );
    }

    // Else if order is greater than degree...
    else if ( order > degree && degree >= 0 )
    {
        // Return zero.
        return 0.0;
    }

    // Else if order and degree are lower than 2...
    else if ( degree <= 1 && order <= 1 )
    {
        // Compute polynomial explicitly.
        return computeGeodesyLegendrePolynomialExplicit( degree, order, polynomialParameter );
    }

    // Else if degree and order are sectoral...
    else if ( degree == order )
    {
        // Obtain polynomial of degree one and order one.
        double degreeOneOrderOnePolynomial = geodesyLegendreCache.getOrElseUpdate(
                    1, 1, polynomialParameter, &computeGeodesyLegendrePolynomial );

        // Obtain prior sectoral polynomial.
        double priorSectoralPolynomial = geodesyLegendreCache.getOrElseUpdate(
                    degree - 1, order - 1, polynomialParameter, &computeGeodesyLegendrePolynomial );

        // Compute polynomial.
        return computeGeodesyLegendrePolynomialDiagonal(
                    degree, degreeOneOrderOnePolynomial, priorSectoralPolynomial );
    }

    // Else degree and order are zonal/tessoral...
    else
    {
        // Obtain prior degree polynomial.
        double oneDegreePriorPolynomial = geodesyLegendreCache.getOrElseUpdate(
                    degree - 1, order, polynomialParameter, &computeGeodesyLegendrePolynomial );

        // Obtain two degrees prior polynomial.
        double twoDegreesPriorPolynomial = geodesyLegendreCache.getOrElseUpdate(
                    degree - 2, order, polynomialParameter, &computeGeodesyLegendrePolynomial );

        // Compute polynomial.
        return computeGeodesyLegendrePolynomialVertical( degree,
                                                         order,
                                                         polynomialParameter,
                                                         oneDegreePriorPolynomial,
                                                         twoDegreesPriorPolynomial );
    }
}

//! Compute derivative of unnormalized Legendre polynomial.
double computeLegendrePolynomialDerivative( const int order,
                                            const double polynomialParameter,
                                            const double currentLegendrePolynomial,
                                            const double incrementedLegendrePolynomial )
{
    // Return polynomial derivative.
    return incrementedLegendrePolynomial
            / std::sqrt( 1.0 - polynomialParameter * polynomialParameter )
            - static_cast< double >( order ) * polynomialParameter
            / ( 1.0 - polynomialParameter * polynomialParameter )
            * currentLegendrePolynomial;
}

//! Compute derivative of geodesy-normalized Legendre polynomial.
double computeGeodesyLegendrePolynomialDerivative( const int degree,
                                                   const int order,
                                                   const double polynomialParameter,
                                                   const double currentLegendrePolynomial,
                                                   const double incrementedLegendrePolynomial )
{
    // Compute normalization correction factor.
    double normalizationCorrection = std::sqrt( ( static_cast< double >( degree )
                                                  + static_cast< double >( order ) + 1.0 )
                                                * ( static_cast< double >( degree - order ) ) );

    // If order is zero apply multiplication factor.
    if ( order == 0 )
    {
        normalizationCorrection *= std::sqrt( 0.5 );
    }

    // Return polynomial derivative.
    return normalizationCorrection * incrementedLegendrePolynomial
            / std::sqrt( 1.0 - polynomialParameter * polynomialParameter )
            - static_cast< double >( order ) * polynomialParameter
            / ( 1.0 - polynomialParameter * polynomialParameter )
            * currentLegendrePolynomial;
}

//! Compute low degree/order unnormalized Legendre polynomial explicitly.
double computeLegendrePolynomialExplicit( const int degree,
                                          const int order,
                                          const double polynomialParameter )
{
    // If 0,0 term is requested return polynomial value.
    if ( degree == 0 && order == 0 )
    {
        return 1.0;
    }

    // Else if 1,0 term is requested return polynomial value.
    else if ( degree == 1 && order == 0 )
    {
        return polynomialParameter;
    }

    // Else if 1,1 term is requested return Legendre polynomial value.
    else if ( degree == 1 && order == 1 )
    {
        return std::sqrt( 1 - polynomialParameter * polynomialParameter );
    }

    // Else the requested term cannot be computed; throw a run-time error.
    else
    {
        // Set error message.
        std::stringstream errorMessage;
        errorMessage << "Error: computation of Legendre polynomial of = " << degree
                     << " and order = " << order << " is not supported." << std::endl;

        // Throw a run-time error.
        boost::throw_exception( boost::enable_error_info( std::runtime_error(
               errorMessage.str( ) ) ) );
    }
}

//! Compute low degree/order geodesy-normalized Legendre polynomials explicitly.
double computeGeodesyLegendrePolynomialExplicit( const int degree,
                                                 const int order,
                                                 const double polynomialParameter )
{
    // If 0,0 term is requested return Legendre polynomial value.
    if ( degree == 0 && order == 0 )
    {
        return 1.0;
    }

    // Else if 1,0 term is requested return polynomial value.
    else if ( degree == 1 && order == 0 )
    {
        return std::sqrt( 3.0 ) * polynomialParameter;
    }

    // Else if 1,1 term is requested return polynomial value.
    else if ( degree == 1 && order == 1 )
    {
        return std::sqrt( 3.0 - 3.0 * polynomialParameter * polynomialParameter );
    }

    // Else the requested term cannot be computed; throw a run-time error.
    else
    {
        // Set error message.
        std::stringstream errorMessage;
        errorMessage << "Error: computation of Legendre polynomial of = " << degree
                     << " and order = " << order << " is not supported." << std::endl;

        // Throw a run-time error.
        boost::throw_exception( boost::enable_error_info( std::runtime_error(
               errorMessage.str( ) ) ) );
    }
}

//! Compute unnormalized Legendre polynomial through sectoral recursion.
double computeLegendrePolynomialDiagonal( const int degree,
                                          const double degreeOneOrderOnePolynomial,
                                          const double priorSectoralPolynomial )
{
    // Return polynomial.
    return ( 2.0 * static_cast< double >( degree ) - 1.0 )
            * degreeOneOrderOnePolynomial * priorSectoralPolynomial;

}

//! Compute geodesy-normalized Legendre polynomial through sectoral recursion.
double computeGeodesyLegendrePolynomialDiagonal( const int degree,
                                                 const double degreeOneOrderOnePolynomial,
                                                 const double priorSectoralPolynomial )
{
    // Return polynomial.
    return std::sqrt( ( 2.0 * static_cast< double >( degree ) + 1.0 )
                      / ( 6.0 * static_cast< double >( degree ) ) )
            * degreeOneOrderOnePolynomial * priorSectoralPolynomial;
}

//! Compute unnormalized Legendre polynomial through degree recursion.
double computeLegendrePolynomialVertical( const int degree,
                                          const int order,
                                          const double polynomialParameter,
                                          const double oneDegreePriorPolynomial,
                                          const double twoDegreesPriorPolynomial )
{
    // Return polynomial.
    return ( ( 2.0 * static_cast< double >( degree ) - 1.0 ) * polynomialParameter
             * oneDegreePriorPolynomial - ( static_cast< double >( degree + order ) - 1.0 )
             * twoDegreesPriorPolynomial ) / ( static_cast< double >( degree - order ) );
}

//! Compute geodesy-normalized Legendre polynomial through degree recursion.
double computeGeodesyLegendrePolynomialVertical( const int degree,
                                                 const int order,
                                                 const double polynomialParameter,
                                                 const double oneDegreePriorPolynomial,
                                                 const double twoDegreesPriorPolynomial )
{
    // Return polynomial.
    return std::sqrt( ( 2.0 * static_cast< double >( degree ) + 1.0 )
                      / ( static_cast< double >( degree + order ) )
                      / ( static_cast< double >( degree - order ) ) )
            * ( std::sqrt( 2.0 * static_cast< double >( degree ) - 1.0 ) * polynomialParameter
                * oneDegreePriorPolynomial
                - std::sqrt( ( static_cast< double >( degree + order ) - 1.0 )
                             * ( static_cast< double >( degree - order ) - 1.0 )
                             / ( 2.0 * static_cast< double >( degree ) - 3.0 ) )
                * twoDegreesPriorPolynomial );
}

//! Define overloaded 'equals' operator for use with 'Point' structure.
bool operator==( const Point& polynomialArguments1, const Point& polynomialArguments2 )
{
    bool equal = polynomialArguments1.degree == polynomialArguments2.degree
            && polynomialArguments1.order == polynomialArguments2.order
            && polynomialArguments1.polynomialParameter
            == polynomialArguments2.polynomialParameter;

    std::cout << "Point "<< writeLegendrePolynomialStructureToString(
                     polynomialArguments1 ).c_str( )
              << " == " << writeLegendrePolynomialStructureToString( polynomialArguments2 ).c_str( )
              << " : " << ( equal ? "TRUE" : "FALSE" ) << std::endl;

    return equal;
}

//! Set hash value.
std::size_t hash_value( Point const& polynomialArguments )
{
    std::size_t seed = 0;
    boost::hash_combine( seed, polynomialArguments.degree );
    boost::hash_combine( seed, polynomialArguments.order );
    boost::hash_combine( seed, polynomialArguments.polynomialParameter );

    std::cout << "Calculated hash for point "<< writeLegendrePolynomialStructureToString(
                     polynomialArguments ).c_str( ) << " as " << seed << std::endl;

    return seed;
}

//! Get Legendre polynomial from cache when possible, and from direct computation otherwise.
double LegendreCache::getOrElseUpdate(
        const int degree, const int order, const double polynomialParameter,
        const LegendrePolynomialFunction legendrePolynomialFunction )
{
    // Initialize structure with polynomial arguments.
    Point polynomialArguments( degree, order, polynomialParameter );

    std::cout << "Query for entry: "
              << writeLegendrePolynomialStructureToString( polynomialArguments ).c_str( )
              << std::endl;

    // Initialize cache iterator.
    CacheTable::iterator cachedEntry = backendCache.find( polynomialArguments );

    // If the requested polynomial was not found in cache, compute polynomial.
    if ( cachedEntry == backendCache.end( ) )
    {
        double legendrePolynomial = legendrePolynomialFunction( degree, order,
                                                                polynomialParameter );
        // If cache is full, remove the oldest element.
        if ( history.full( ) )
        {
            std::cout << "History is full, removing oldest element: "
                      << writeLegendrePolynomialStructureToString( history[ 0 ] ).c_str( )
                      << std::endl;
            std::cout << "BEFORE:" << std::endl;
            dumpLegendrePolynomialCacheData( std::cout, backendCache, history );
            backendCache.erase( backendCache.find( history[ 0 ] ) );
            history.pop_front( );
            std::cout << "AFTER:" << std::endl;
            dumpLegendrePolynomialCacheData( std::cout, backendCache, history );
        }

        std::cout << "Inserting new point in history: "
                  << writeLegendrePolynomialStructureToString( polynomialArguments ).c_str( )
                  << " => " << legendrePolynomial << std::endl;
        std::cout << "BEFORE:" << std::endl;
        dumpLegendrePolynomialCacheData( std::cout, backendCache, history );

        // Insert computed polynomial into cache.
        backendCache.insert( std::pair< Point, double >( polynomialArguments,
                                                         legendrePolynomial ) );
        history.push_back( polynomialArguments );

        std::cout << "AFTER:" << std::endl;
        dumpLegendrePolynomialCacheData( std::cout, backendCache, history );

        // Return polynomial value from computation.
        return legendrePolynomial;
    }

    // Else the requested polynomial was found in cache; return polynomial value from cache entry.
    else
    {
        std::cout << "Found entry in cache: "
                  << writeLegendrePolynomialStructureToString( cachedEntry->first ).c_str( )
                  << " => "
                  << cachedEntry->second
                  << std::endl;
        return cachedEntry->second;
    }
}

//! Initialize LegendreCache objects.
LegendreCache::LegendreCache( ) : history( MAXIMUM_CACHE_ENTRIES ) { }

//! Write contents of Legendre polynomial structure to string.
std::string writeLegendrePolynomialStructureToString( const Point legendrePolynomialStructure )
{
    std::stringstream buffer;

    buffer << "(" << legendrePolynomialStructure.degree << ", "
           << legendrePolynomialStructure.order << ", "
           << legendrePolynomialStructure.polynomialParameter << ")";

    return buffer.str( );
}

//! Dump Legendre polynomial cache data to stream (table and history).
void dumpLegendrePolynomialCacheData( std::ostream& outputStream,
                                      boost::unordered_map< Point, double > cacheTable,
                                      boost::circular_buffer< Point > cacheHistory )
{
    outputStream << "Table:\n";

    for ( boost::unordered_map< Point, double >::iterator iteratorCacheTable = cacheTable.begin( );
          iteratorCacheTable != cacheTable.end( ); iteratorCacheTable++ )
    {
        outputStream << "\t" << writeLegendrePolynomialStructureToString(
                            iteratorCacheTable->first ).c_str( ) << " => "
                     << iteratorCacheTable->second << std::endl;
    }

    outputStream << "History:\n";

    for ( boost::circular_buffer< Point >::iterator iteratorCacheHistory = cacheHistory.begin( );
          iteratorCacheHistory != cacheHistory.end( ); iteratorCacheHistory++ )
    {
        outputStream << "\t"
                     << writeLegendrePolynomialStructureToString( *iteratorCacheHistory ).c_str( )
                     << ", ";
    }

    outputStream << std::endl;
}

} // namespace basic_mathematics
} // namespace tudat
