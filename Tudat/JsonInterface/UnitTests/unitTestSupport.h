/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include "Tudat/JsonInterface/Support/deserialization.h"
#include "Tudat/JsonInterface/Support/utilities.h"
#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"

namespace tudat
{

namespace json_interface
{

boost::filesystem::path currentDirectory( )
{
    return boost::filesystem::path( __FILE__ ).parent_path( );
}

boost::filesystem::path inputDirectory( )
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


void checkCloseIntegrationResults( const std::map< double, Eigen::VectorXd >& results1,
                                   const std::map< double, Eigen::VectorXd >& results2,
                                   const std::vector< unsigned int > indices,
                                   const std::vector< unsigned int > sizes,
                                   const double tolerance )
{
    // Check size of maps
    BOOST_CHECK_EQUAL( results1.size( ), results2.size( ) );

    // Check initial epochs
    const double initialEpoch1 = results1.begin( )->first;
    const double initialEpoch2 = results2.begin( )->first;
    BOOST_CHECK_SMALL( std::fabs( initialEpoch1 - initialEpoch2 ), tolerance );

    // Check final epochs
    const double finalEpoch1 = ( --results1.end( ) )->first;
    const double finalEpoch2 = ( --results2.end( ) )->first;
    BOOST_CHECK_SMALL( std::fabs( finalEpoch1 - finalEpoch2 ), tolerance );

    // Check norm of requested vector segments
    const Eigen::VectorXd initialState1 = results1.begin( )->second;
    const Eigen::VectorXd initialState2 = results2.begin( )->second;
    const Eigen::VectorXd finalState1 = ( --results1.end( ) )->second;
    const Eigen::VectorXd finalState2 = ( --results2.end( ) )->second;
    for ( unsigned int i = 0; i < indices.size( ); ++i )
    {
        // Initial step
        const double initialNorm1 = initialState1.segment( indices.at( i ), sizes.at( i ) ).norm( );
        const double initialNorm2 = initialState2.segment( indices.at( i ), sizes.at( i ) ).norm( );
        BOOST_CHECK_CLOSE_FRACTION( initialNorm1, initialNorm2, tolerance );

        // Final step
        const double finalNorm1 = finalState1.segment( indices.at( i ), sizes.at( i ) ).norm( );
        const double finalNorm2 = finalState2.segment( indices.at( i ), sizes.at( i ) ).norm( );
        BOOST_CHECK_CLOSE_FRACTION( finalNorm1, finalNorm2, tolerance );
    }
}

#define BOOST_CHECK_CLOSE_INTEGRATION_RESULTS( results1, results2, indices, sizes, tolerance ) \
    tudat::json_interface::checkCloseIntegrationResults( results1, results2, indices, sizes, tolerance )

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_UNITTESTSUPPORT
