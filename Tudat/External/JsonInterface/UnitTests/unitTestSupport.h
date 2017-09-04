/*    Copyright (c) 2010-2017, Delft University of Technology
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

#include "Tudat/External/JsonInterface/Support/modular.h"
#include "Tudat/External/JsonInterface/Support/utilities.h"

#include <Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h>

namespace tudat
{

namespace json_interface
{

path currentDirectory( )
{
    return path( __FILE__ ).parent_path( );
}

path inputDirectory( )
{
    return currentDirectory( ) / "inputs";
}

path outputDirectory( )
{
    return currentDirectory( ) / "outputs";
}

template< typename T = json >
T readInputFile( const std::string& filename, const std::string& extension = "json" )
{
    const path filePath = inputDirectory( ) / ( filename + "." + extension );
    boost::filesystem::current_path( filePath.parent_path( ) );
    return readJSON( filePath.string( ) ).get< T >( );
}


template< typename T >
void checkJsonEquivalent( const T& left, const T& right )
{
    const std::string fromFile = json( left ).dump( 2 );
    const std::string manual = json( right ).dump( 2 );
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
    const std::vector< Enum > values = readInputFile< std::vector< Enum > >( filename );
    BOOST_CHECK_EQUAL_JSON( values, supportedValues );
}

#define BOOST_CHECK_EQUAL_ENUM( filename, stringValues, usupportedValues ) \
    tudat::json_interface::checkConsistentEnum( filename, stringValues, usupportedValues )


void checkCloseIntegrationResults( const std::map< double, Eigen::VectorXd >& results1,
                                   const std::map< double, Eigen::VectorXd >& results2,
                                   const double tolerance )
{
    // Results 1

    const double initialEpoch1 = results1.begin( )->first;
    const Eigen::VectorXd initialState1 = results1.begin( )->second;

    const double finalEpoch1 = ( --results1.end( ) )->first;
    const Eigen::VectorXd finalState1 = ( --results1.end( ) )->second;


    // Results 2

    const double initialEpoch2 = results2.begin( )->first;
    const Eigen::VectorXd initialState2 = results2.begin( )->second;

    const double finalEpoch2 = ( --results2.end( ) )->first;
    const Eigen::VectorXd finalState2 = ( --results2.end( ) )->second;


    // Compare initial conditions

    BOOST_CHECK_SMALL( std::fabs( initialEpoch1 - initialEpoch2 ), tolerance );
    BOOST_CHECK_SMALL( ( initialState1 - initialState2 ).norm( ), tolerance );


    // Compare final conditions

    BOOST_CHECK_SMALL( std::fabs( finalEpoch1 - finalEpoch2 ), tolerance );
    BOOST_CHECK_SMALL( ( finalState1 - finalState2 ).norm( ), tolerance );
}

#define BOOST_CHECK_CLOSE_INTEGRATION_RESULTS( results1, results2, tolerance ) \
    tudat::json_interface::checkCloseIntegrationResults( results1, results2, tolerance )

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_UNITTESTSUPPORT
