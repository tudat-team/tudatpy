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


void checkCloseIntegrationResults(
        const boost::shared_ptr< propagators::SingleArcDynamicsSimulator< double, double > >& jsonSimulator,
        const boost::shared_ptr< propagators::SingleArcDynamicsSimulator< double, double > >& manualSimulator,
        const double tolerance )
{
    // JSON

    const std::map< double, Eigen::VectorXd > integrationResultJSON =
            jsonSimulator->getEquationsOfMotionNumericalSolution( );

    const double initialEpochJSON = integrationResultJSON.begin( )->first;
    const Eigen::VectorXd initialStateJSON = integrationResultJSON.begin( )->second;

    const double finalEpochJSON = ( --integrationResultJSON.end( ) )->first;
    const Eigen::VectorXd finalStateJSON = ( --integrationResultJSON.end( ) )->second;


    // Manual

    const std::map< double, Eigen::VectorXd > integrationResult =
            manualSimulator->getEquationsOfMotionNumericalSolution( );

    const double initialEpoch = integrationResult.begin( )->first;
    const Eigen::VectorXd initialState = integrationResult.begin( )->second;

    const double finalEpoch = ( --integrationResult.end( ) )->first;
    const Eigen::VectorXd finalState = ( --integrationResult.end( ) )->second;


    // Compare initial conditions

    BOOST_CHECK_SMALL( std::fabs( initialEpochJSON - initialEpoch ), tolerance );
    BOOST_CHECK_SMALL( ( initialStateJSON - initialState ).norm( ), tolerance );


    // Compare final conditions

    BOOST_CHECK_SMALL( std::fabs( finalEpochJSON - finalEpoch ), tolerance );
    BOOST_CHECK_SMALL( ( finalStateJSON - finalState ).norm( ), tolerance );
}

#define BOOST_CHECK_CLOSE_INTEGRATION_RESULTS( jsonSimulator, manualSimulator, tolerance ) \
    tudat::json_interface::checkCloseIntegrationResults( jsonSimulator, manualSimulator, tolerance )

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_UNITTESTSUPPORT
