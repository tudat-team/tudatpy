/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_UNITTESTSUPPORT
#define TUDAT_JSONINTERFACE_UNITTESTSUPPORT

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/interface/json/support/deserialization.h"
#include "tudat/interface/json/support/utilities.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"

namespace tudat
{

namespace json_interface
{

inline boost::filesystem::path currentDirectory( )
{
    return boost::filesystem::path( __FILE__ ).parent_path( );
}

inline boost::filesystem::path inputDirectory( )
{
    boost::filesystem::path matlabInputDirectory = currentDirectory( ) / "matlab_inputs";
    if ( boost::filesystem::exists( matlabInputDirectory ) )
    {
        return matlabInputDirectory;
    }
    else
    {
        return currentDirectory( ) / "inputs";
    }
}

template< typename T = nlohmann::json >
T parseJSONFile( std::string file )
{
    if ( boost::filesystem::path( file ).extension( ).empty( ) )
    {
        file += ".json";
    }
    boost::filesystem::current_path( boost::filesystem::path( file ).parent_path( ) );
    return readJSON( file );
}


template< typename T >
void checkJsonEquivalent( const T& left, const T& right )
{
    const std::string fromFile = nlohmann::json( left ).dump( 2 );
    const std::string manual = nlohmann::json( right ).dump( 2 );
    BOOST_CHECK_EQUAL( fromFile, manual );
}

#define BOOST_CHECK_EQUAL_JSON( left, right ) tudat::json_interface::checkJsonEquivalent( left, right )


template< typename Enum >
void checkConsistentEnum( const std::string& filename,
                          const std::map< Enum, std::string >& stringValues,
                          const std::vector< Enum >& usupportedValues )
{
    using namespace json_interface;

    // Create vector of supported values
    std::vector< Enum > supportedValues;
    for ( auto entry : stringValues )
    {
        Enum value = entry.first;
        if ( ! contains( usupportedValues, value ) )
        {
            supportedValues.push_back( value );
        }
    }

    // Check that values and supportedValues are equivalent
    const std::vector< Enum > values = parseJSONFile< std::vector< Enum > >( filename );
    BOOST_CHECK( containsAllOf( values, supportedValues ) && containsAllOf( supportedValues, values ) );
}

#define BOOST_CHECK_EQUAL_ENUM( filename, stringValues, usupportedValues ) \
    tudat::json_interface::checkConsistentEnum( filename, stringValues, usupportedValues )

template< typename ContainerType >
void checkCloseIntegrationResultsMatrix( const std::map< double, ContainerType >& results1,
                                         const std::map< double, ContainerType >& results2,
                                         const std::vector< std::pair< unsigned int, unsigned int > > indices,
                                         const std::vector< std::pair< unsigned int, unsigned int > > sizes,
                                         const std::vector< double > absoluteTolerances )
{
    // Check size of maps
    BOOST_CHECK_EQUAL( results1.size( ), results2.size( ) );

    auto it1 = results1.begin( );
    auto it2 = results2.begin( );

    for( unsigned int i = 0; i < results1.size( ); i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( it1->first - it2->first ),
                           2.0 * it1->first * std::numeric_limits< double >::epsilon( ) );
        const Eigen::MatrixXd state1 = it1->second;
        const Eigen::MatrixXd state2 = it2->second;

        for ( unsigned int i = 0; i < indices.size( ); ++i )
        {
            for( unsigned int row = 0; row < sizes.at( i ).first; row++ )
            {
                for( unsigned int col = 0; col < sizes.at( i ).second; col++ )
                {
                    BOOST_CHECK_SMALL(
                                std::fabs(
                                    state1( indices.at( i ).first + row, indices.at( i ).second + col ) -
                                    state2( indices.at( i ).first + row, indices.at( i ).second + col ) ), absoluteTolerances.at( i ) );
                }
            }
        }

    }
}

inline void checkCloseIntegrationResults( const std::map< double, Eigen::VectorXd >& results1,
                                   const std::map< double, Eigen::VectorXd >& results2,
                                   const std::vector< unsigned int > indices,
                                   const std::vector< unsigned int > sizes,
                                   const double tolerance )
{
    std::vector< std::pair< unsigned int, unsigned int > > indicesFull;
    std::vector< std::pair< unsigned int, unsigned int > > sizesFull;
    std::vector< double > absoluteTolerancesFull;

    for( unsigned int i = 0 ; i < indices.size( ); i++ )
    {
        indicesFull.push_back( std::make_pair( indices.at( i ), 0 ) ) ;
        sizesFull.push_back( std::make_pair( sizes.at( i ), 1 ) ) ;
        absoluteTolerancesFull.push_back( tolerance );
    }
    checkCloseIntegrationResultsMatrix< Eigen::VectorXd >( results1, results2, indicesFull, sizesFull, absoluteTolerancesFull );
}

#define BOOST_CHECK_CLOSE_INTEGRATION_RESULTS( results1, results2, indices, sizes, tolerance ) \
    tudat::json_interface::checkCloseIntegrationResults( results1, results2, indices, sizes, tolerance )

#define BOOST_CHECK_CLOSE_INTEGRATION_MATRIX_RESULTS( results1, results2, indices, sizes, tolerance ) \
    tudat::json_interface::checkCloseIntegrationResultsMatrix< Eigen::MatrixXd >( results1, results2, indices, sizes, tolerance );

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_UNITTESTSUPPORT
