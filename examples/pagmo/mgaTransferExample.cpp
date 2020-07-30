/*    Copyright (c) 2010-2019, Delft University of Technology
*    All rigths reserved
*
*    This file is part of the Tudat. Redistribution and use in source and
*    binary forms, with or without modification, are permitted exclusively
*    under the terms of the Modified BSD license. You should have received
*    a copy of the license with this file. If not, please or visit:
*    http://tudat.tudelft.nl/LICENSE.
*/

#include <iostream>
#include <fstream>

#include <boost/filesystem.hpp>

#include "pagmo/island.hpp"
#include "pagmo/io.hpp"
#include "pagmo/problem.hpp"
#include <pagmo/rng.hpp>
#include "pagmo/algorithms/sade.hpp"
#include "problems/multipleGravityAssist.h"
#include "tudat/astro/ephemerides/approximatePlanetPositions.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/astro/mission_segments/createTransferTrajectory.h"
#include "tudat/simulation/environment/body.h"

using namespace pagmo;
using namespace tudat;
using namespace ephemerides;
using namespace gravitation;
using namespace simulation_setup;

simulation_setup::NamedBodyMap getApproximatePlanetBodyMap( )
{


    NamedBodyMap bodyMap;
    bodyMap[ "Sun" ] = std::make_shared< Body >( );
    bodyMap[ "Mercury" ] = std::make_shared< Body >( );
    bodyMap[ "Venus" ] = std::make_shared< Body >( );
    bodyMap[ "Earth" ] = std::make_shared< Body >( );
    bodyMap[ "Mars" ] = std::make_shared< Body >( );
    bodyMap[ "Jupiter" ] = std::make_shared< Body >( );
    bodyMap[ "Saturn" ] = std::make_shared< Body >( );

    bodyMap[ "Sun" ]->setEphemeris( std::make_shared< ConstantEphemeris >( Eigen::Vector6d::Zero( ) ) );
    bodyMap[ "Mercury" ]->setEphemeris( std::make_shared< ApproximatePlanetPositions >(
                                            ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mercury ) );
    bodyMap[ "Venus" ]->setEphemeris( std::make_shared< ApproximatePlanetPositions >(
                                          ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus ) );
    bodyMap[ "Earth" ]->setEphemeris( std::make_shared< ApproximatePlanetPositions >(
                                          ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter ) );
    bodyMap[ "Mars" ]->setEphemeris( std::make_shared< ApproximatePlanetPositions >(
                                          ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter ) );
    bodyMap[ "Jupiter" ]->setEphemeris( std::make_shared< ApproximatePlanetPositions >(
                                            ApproximatePlanetPositionsBase::BodiesWithEphemerisData::jupiter ) );
    bodyMap[ "Saturn" ]->setEphemeris( std::make_shared< ApproximatePlanetPositions >(
                                           ApproximatePlanetPositionsBase::BodiesWithEphemerisData::saturn ) );

    bodyMap[ "Sun" ]->setGravityFieldModel( std::make_shared< GravityFieldModel >( 1.32712428e20 ) );
    bodyMap[ "Mercury" ]->setGravityFieldModel( std::make_shared< GravityFieldModel >( 2.2321e13 ) );
    bodyMap[ "Venus" ]->setGravityFieldModel( std::make_shared< GravityFieldModel >( 3.24860e14 ) );
    bodyMap[ "Earth" ]->setGravityFieldModel( std::make_shared< GravityFieldModel >( 3.9860119e14 ) );
    bodyMap[ "Mars" ]->setGravityFieldModel( std::make_shared< GravityFieldModel >( 4.282837e13 ) );
    bodyMap[ "Jupiter" ]->setGravityFieldModel( std::make_shared< GravityFieldModel >( 1.267e17 ) );
    bodyMap[ "Saturn" ]->setGravityFieldModel( std::make_shared< GravityFieldModel >( 3.79e16 ) );

    return bodyMap;
}

namespace std
{
template <typename T>
ostream& operator<<(ostream& output, std::vector<T> const& values)
{
    std::cout<<"[ ";
    for (auto const& value : values)
    {
        output << value <<", ";
    }
    std::cout<<" ]"<<std::endl;
    return output;
}
}
//! Execute  main
int main( )
{
    //Set seed for reproducible results
    pagmo::random_device::set_seed( 123456789 );

    // Set transfer order
    std::vector< std::string > bodyOrder = {
        "Earth", "Venus", "Venus", "Earth", "Jupiter", "Saturn" };
    int numberOfNodes = bodyOrder.size( );

    // Create leg settings (all unpowered)
    std::vector< std::shared_ptr< TransferLegSettings > > transferLegSettings;
    transferLegSettings.resize( numberOfNodes - 1 );
    transferLegSettings[ 0 ] = unpoweredLeg( );
    transferLegSettings[ 1 ] = dsmVelocityBasedLeg( );
    transferLegSettings[ 2 ] = dsmVelocityBasedLeg( );
    transferLegSettings[ 3 ] = unpoweredLeg( );
    transferLegSettings[ 4 ] = unpoweredLeg( );

    // Define minimum periapsis altitudes for flybys;
    std::map< std::string, double >minimumPeriapses;
    minimumPeriapses[ "Venus" ] = 6351800.0;
    minimumPeriapses[ "Earth" ] = 6778000.0;
    minimumPeriapses[ "Mars" ] = 3696200.0;
    minimumPeriapses[ "Jupiter" ] =  600000000;

    std::vector< std::shared_ptr< TransferNodeSettings > > transferNodeSettings;
    transferNodeSettings.resize( numberOfNodes );
    transferNodeSettings[ 0 ] = escapeAndCaptureNode( std::numeric_limits< double >::infinity( ), 0.0 );
    transferNodeSettings[ 1 ] = swingbyNode( minimumPeriapses.at( bodyOrder.at( 1 ) ) );
    transferNodeSettings[ 2 ] = swingbyNode( minimumPeriapses.at( bodyOrder.at( 2 ) ) );
    transferNodeSettings[ 3 ] = swingbyNode( minimumPeriapses.at( bodyOrder.at( 3 ) ) );
    transferNodeSettings[ 4 ] = swingbyNode( minimumPeriapses.at( bodyOrder.at( 4 ) ) );
    transferNodeSettings[ 5 ] = captureAndInsertionNode( 1.0895e8 / 0.02, 0.98 );

    printTransferParameterDefinition( transferLegSettings, transferNodeSettings );
    simulation_setup::NamedBodyMap bodyMap = getApproximatePlanetBodyMap( );

    // Define search bounds: first parameter is start date, following parameters are leg durations
    std::vector< std::vector< double > > bounds( 2, std::vector< double >( 14, 0.0 ) );
    bounds[ 0 ][ 0 ] = -2000.0; //MJD2000
    bounds[ 1 ][ 0 ] = 0.0; //MJD2000
    bounds[ 0 ][ 1 ] = 50.0;
    bounds[ 1 ][ 1 ] = 500.0;
    bounds[ 0 ][ 2 ] = 100.0;
    bounds[ 1 ][ 2 ] = 500.0;
    bounds[ 0 ][ 3 ] = 50.0;
    bounds[ 1 ][ 3 ] = 500.0;
    bounds[ 0 ][ 4 ] = 500.0;
    bounds[ 1 ][ 4 ] = 2000.0;
    bounds[ 0 ][ 5 ] = 1000.0;
    bounds[ 1 ][ 5 ] = 10000.0;
//    bounds[ 0 ][ 6 ] = 0.0;
//    bounds[ 1 ][ 6 ] = 2.0E3;;
//    bounds[ 0 ][ 7 ] = 0.0;
//    bounds[ 1 ][ 7 ] = 2.0 * mathematical_constants::PI;

//    bounds[ 0 ][ 8 ] = 0.0;
//    bounds[ 1 ][ 8 ] = 2.0 * mathematical_constants::PI;

//    bounds[ 0 ][ 9 ] = 0.05;
//    bounds[ 1 ][ 9 ] = 0.95;
    int currentIndex = 5;
    for( int i = 1; i < 3; i++ )
    {
        // periapsis bounds
        bounds[ 0 ][ currentIndex + ( i - 1 ) * 4 + 1 ] = minimumPeriapses.at( bodyOrder.at( i ) );
        bounds[ 1 ][ currentIndex + ( i - 1 ) * 4 + 1 ] = 1.5 * minimumPeriapses.at( bodyOrder.at( i ) );

        // orbit orientation bounds
        bounds[ 0 ][ currentIndex + ( i - 1 ) * 4 + 2 ] = 0.0;
        bounds[ 1 ][ currentIndex + ( i - 1 ) * 4 + 2 ] = 2.0 * mathematical_constants::PI;

        // Swingby Delta V bounds
        bounds[ 0 ][ currentIndex + ( i - 1 ) * 4 + 3 ] = 0;
        bounds[ 1 ][ currentIndex + ( i - 1 ) * 4 + 3 ] = 2500.0;

        // DSM TOF fraction
        bounds[ 0 ][ currentIndex + ( i - 1 ) * 4 + 4 ] = 0.05;
        bounds[ 1 ][ currentIndex + ( i - 1 ) * 4 + 4 ] = 0.95;
    }

    // Create object to compute the problem fitness
    problem prob{ MultipleGravityAssist(
                    bodyMap, transferLegSettings, transferNodeSettings, bodyOrder, "Sun", bounds ) };


    // Select NSGA2 algorithm for priblem
    algorithm algo{sade( )};

    // Create an island with 1000 individuals
    island isl{algo, prob, 500 };

    // Evolve for 512 generations
    for( int i = 0 ; i < 10000; i++ )
    {

        isl.evolve( );
        while( isl.status( ) != pagmo::evolve_status::idle &&
               isl.status( ) != pagmo::evolve_status::idle_error )
        {
            isl.wait( );
        }

        isl.wait_check( ); // Raises errors

        if( i% 100 == 0 )
        {
            std::cout<<i<<" "<<isl.get_population().champion_f()[0]<<"***  "<<
                       isl.get_population().champion_x()<<std::endl;        // Write current iteration results to file
        }
//        printPopulationToFile( isl.get_population( ).get_x( ), "mo_mga_EVEEJ_" + std::to_string( i ), false );
//        printPopulationToFile( isl.get_population( ).get_f( ), "mo_mga_EVEEJ_" + std::to_string( i ), true );
    }

    return 0;

}
