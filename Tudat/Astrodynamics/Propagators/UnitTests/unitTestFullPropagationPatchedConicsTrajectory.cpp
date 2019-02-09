/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <iomanip>

#include <Tudat/SimulationSetup/tudatEstimationHeader.h>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <Eigen/Core>
#include "Tudat/Basics/testMacros.h"

#include "Tudat/SimulationSetup/PropagationSetup/propagationLambertTargeterFullProblem.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationPatchedConicFullProblem.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/trajectory.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/exportTrajectory.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/planetTrajectory.h"


namespace tudat
{

namespace unit_tests
{

using namespace tudat;

//! Test of the full propagation of a trajectory
BOOST_AUTO_TEST_SUITE( testFullPropagationTrajectory )

//! Test of the full propagation of a MGA trajectory
BOOST_AUTO_TEST_CASE( testFullPropagationMGA )
{

    std::cout.precision(20);

    //    std::cout << "Cassini trajectory: " << "\n\n";


    // Specify number and type of legs.
    int numberOfLegs = 6;
    std::vector< transfer_trajectories::TransferLegType > legTypeVector;
    legTypeVector.resize( numberOfLegs );
    legTypeVector[ 0 ] = transfer_trajectories::mga_Departure;
    legTypeVector[ 1 ] = transfer_trajectories::mga_Swingby;
    legTypeVector[ 2 ] = transfer_trajectories::mga_Swingby;
    legTypeVector[ 3 ] = transfer_trajectories::mga_Swingby;
    legTypeVector[ 4 ] = transfer_trajectories::mga_Swingby;
    legTypeVector[ 5 ] = transfer_trajectories::capture;

    // Name of the bodies involved in the trajectory
    std::vector< std::string > nameBodiesTrajectory;
    nameBodiesTrajectory.push_back("Earth");
    nameBodiesTrajectory.push_back("Venus");
    nameBodiesTrajectory.push_back("Mars");
    nameBodiesTrajectory.push_back("Earth");
    nameBodiesTrajectory.push_back("Jupiter");
    nameBodiesTrajectory.push_back("Saturn");



    std::vector< std::string > centralBody; centralBody.push_back( "Sun" );
    std::string bodyToPropagate = "spacecraft";


    spice_interface::loadStandardSpiceKernels( );


    // Define gravitational parameter for each transfer body.
    std::vector< double > gravitationalParametersTransferBodies;
    gravitationalParametersTransferBodies.push_back( simulation_setup::createGravityFieldModel( simulation_setup::getDefaultGravityFieldSettings(
                                                                                                    "Earth", TUDAT_NAN, TUDAT_NAN), "Earth")->getGravitationalParameter());
    gravitationalParametersTransferBodies.push_back( simulation_setup::createGravityFieldModel( simulation_setup::getDefaultGravityFieldSettings(
                                                                                                    "Venus", TUDAT_NAN, TUDAT_NAN), "Venus")->getGravitationalParameter());
    gravitationalParametersTransferBodies.push_back( simulation_setup::createGravityFieldModel( simulation_setup::getDefaultGravityFieldSettings(
                                                                                                    "Venus", TUDAT_NAN, TUDAT_NAN), "Venus")->getGravitationalParameter());
    gravitationalParametersTransferBodies.push_back( simulation_setup::createGravityFieldModel( simulation_setup::getDefaultGravityFieldSettings(
                                                                                                    "Earth", TUDAT_NAN, TUDAT_NAN), "Earth")->getGravitationalParameter());
    gravitationalParametersTransferBodies.push_back( simulation_setup::createGravityFieldModel( simulation_setup::getDefaultGravityFieldSettings(
                                                                                                    "Jupiter", TUDAT_NAN, TUDAT_NAN), "Jupiter")->getGravitationalParameter());
    gravitationalParametersTransferBodies.push_back( simulation_setup::createGravityFieldModel( simulation_setup::getDefaultGravityFieldSettings(
                                                                                                    "Saturn", TUDAT_NAN, TUDAT_NAN), "Saturn")->getGravitationalParameter());


    // Create variable vector.
    std::vector< double > variableVector;
    variableVector.push_back(-789.8117 * physical_constants::JULIAN_DAY); variableVector.push_back(158.302027105278 * physical_constants::JULIAN_DAY);
    variableVector.push_back(449.385873819743 * physical_constants::JULIAN_DAY); variableVector.push_back(54.7489684339665 * physical_constants::JULIAN_DAY);
    variableVector.push_back(1024.36205846918 * physical_constants::JULIAN_DAY); variableVector.push_back(4552.30796805542 * physical_constants::JULIAN_DAY);
    variableVector.push_back(1.0 * physical_constants::JULIAN_DAY);

    double initialTime = variableVector[0];



    // Define ephemerides for the transfer bodies
    std::vector< ephemerides::EphemerisPointer > ephemerisVectorTransferBodies;
    ephemerides::EphemerisPointer ephemerisEarth = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter) ;
    ephemerides::EphemerisPointer ephemerisVenus = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus) ;
    ephemerides::EphemerisPointer ephemerisJupiter = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::jupiter) ;
    ephemerides::EphemerisPointer ephemerisSaturn = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::saturn) ;
    ephemerisVectorTransferBodies.push_back( ephemerisEarth );
    ephemerisVectorTransferBodies.push_back( ephemerisVenus );
    ephemerisVectorTransferBodies.push_back( ephemerisVenus );
    ephemerisVectorTransferBodies.push_back( ephemerisEarth );
    ephemerisVectorTransferBodies.push_back( ephemerisJupiter );
    ephemerisVectorTransferBodies.push_back( ephemerisSaturn );


    // Create body map.
    simulation_setup::NamedBodyMap bodyMap = propagators::setupBodyMapFromUserDefinedEphemeridesForPatchedConicsTrajectory(centralBody[0],
            bodyToPropagate, nameBodiesTrajectory, ephemerisVectorTransferBodies, gravitationalParametersTransferBodies, "ECLIPJ2000" );

    // Create acceleration map.
    std::vector< basic_astrodynamics::AccelerationMap > accelerationMap = propagators::setupAccelerationMapPatchedConicsTrajectory(
                nameBodiesTrajectory.size(), centralBody[0], bodyToPropagate, bodyMap);

    // Create departure and capture variables.
    std::vector< double > semiMajorAxes;
    semiMajorAxes.push_back( std::numeric_limits< double >::infinity( ) ); semiMajorAxes.push_back( 1.0895e8 / 0.02 );
    std::vector< double > eccentricities;
    eccentricities.push_back( 0.0 ); eccentricities.push_back( 0.98 );

    // Create minimum pericenter radii vector
    std::vector< double > minimumPericenterRadii;
    minimumPericenterRadii.push_back( 6778000.0 ); minimumPericenterRadii.push_back( 6351800.0 ); minimumPericenterRadii.push_back( 6351800.0 );
    minimumPericenterRadii.push_back( 6778000.0 ); minimumPericenterRadii.push_back( 600000000.0 ); minimumPericenterRadii.push_back( 600000000.0 );


    // Define integrator settings.
    double fixedStepSize = 3600.0;
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared < numerical_integrators::IntegratorSettings < > > ( numerical_integrators::rungeKutta4, initialTime, fixedStepSize);


    // Compute difference between patched conics trajectory and full problem at departure and at arrival for each leg.
    for( int terminationType = 0; terminationType < 2; terminationType++ )
    {
        std::map< int, std::map< double, Eigen::Vector6d > > patchedConicsResultForEachLeg;
        std::map< int, std::map< double, Eigen::Vector6d > > fullProblemResultForEachLeg;
        std::map< int, std::map< double, Eigen::VectorXd > > dependentVariableResultForEachLeg;

        propagators::fullPropagationPatchedConicsTrajectory(
                    bodyMap, accelerationMap, nameBodiesTrajectory,
                    centralBody[0], bodyToPropagate, legTypeVector, variableVector, minimumPericenterRadii,
                semiMajorAxes, eccentricities, integratorSettings, patchedConicsResultForEachLeg, fullProblemResultForEachLeg,
                dependentVariableResultForEachLeg, static_cast< bool >( terminationType ) );

        for( auto itr : patchedConicsResultForEachLeg )
        {
            for( auto innerItr : itr.second )
            {
                Eigen::Vector6d stateDifference = patchedConicsResultForEachLeg[ itr.first ][ innerItr.first ] -
                        fullProblemResultForEachLeg[ itr.first ][ innerItr.first ];
                for( int i = 0; i < 3; i++ )
                {
                    BOOST_CHECK_SMALL( std::fabs( stateDifference( i ) ), 1.0 );
                    BOOST_CHECK_SMALL( std::fabs( stateDifference( i + 3 ) ), 1.0E-6 );
                    BOOST_CHECK_SMALL( std::fabs( stateDifference( i ) ), 1.0 );
                    BOOST_CHECK_SMALL( std::fabs( stateDifference( i + 3 ) ), 1.0E-6 );
                }
            }
        }

        int numberOfLegs = patchedConicsResultForEachLeg.size( );
        auto legIterator = patchedConicsResultForEachLeg.begin( );
        auto nextLegIterator = patchedConicsResultForEachLeg.begin( );
        nextLegIterator++;

        for( int i = 0; i < numberOfLegs - 1; i++ )
        {
            double legTimeDifference =
                    nextLegIterator->second.begin( )->first - legIterator->second.rbegin( )->first;
            Eigen::Vector6d legStateDifference =
                    nextLegIterator->second.begin( )->second - legIterator->second.rbegin( )->second;

            if( terminationType == 0 )
            {
                BOOST_CHECK_SMALL( std::fabs( legTimeDifference ), 1.0E-6 );
                if( i != 1 )
                {
                    BOOST_CHECK_SMALL( std::fabs( legStateDifference.segment( 0, 3 ).norm( ) ), 1.0 );
                }
            }
            else
            {
                BOOST_CHECK_EQUAL( std::fabs( legTimeDifference ) > 1.0E5, true );
            }

//            std::cout<<"Time diff. "<<legTimeDifference<<std::endl;

//            std::cout<<std::setprecision( 16 )<<"State. "<<( legIterator->second.rbegin( )->second ).transpose( )<<std::endl; //Wrong
//            std::cout<<std::setprecision( 16 )<<"State. "<<( nextLegIterator->second.begin( )->second ).transpose( )<<std::endl; //Correct
//            std::cout<<std::setprecision( 16 )<<"Venus state. "<<ephemerisVenus->getCartesianState( nextLegIterator->second.begin( )->first ).transpose( )<<std::endl<<std::endl;

//            std::cout<<"State diff. "<<legStateDifference.transpose( )<<std::endl<<std::endl;

            legIterator++;
            nextLegIterator++;
        }
//        std::cout<<std::endl<<std::endl;
    }
}


//! Test of the full propagation of a MGA trajectory including deep-space manoeuvres
BOOST_AUTO_TEST_CASE( testFullPropagationMGAwithDSM )
{

    //    std::cout << "Messenger trajectory: " << "\n\n";

    // Specify number and type of legs.
    int numberOfLegs = 5;
    std::vector< transfer_trajectories::TransferLegType > legTypeVector;
    legTypeVector.resize( numberOfLegs );
    legTypeVector[ 0 ] = transfer_trajectories::mga1DsmVelocity_Departure;
    legTypeVector[ 1 ] = transfer_trajectories::mga1DsmVelocity_Swingby;
    legTypeVector[ 2 ] = transfer_trajectories::mga1DsmVelocity_Swingby;
    legTypeVector[ 3 ] = transfer_trajectories::mga1DsmVelocity_Swingby;
    legTypeVector[ 4 ] = transfer_trajectories::capture;


    // Name of the bodies involved in the trajectory
    std::vector< std::string > nameBodiesAndManoeuvresTrajectory;
    nameBodiesAndManoeuvresTrajectory.push_back("Earth");
    nameBodiesAndManoeuvresTrajectory.push_back("DSM1");
    nameBodiesAndManoeuvresTrajectory.push_back("Earth");
    nameBodiesAndManoeuvresTrajectory.push_back("DSM2");
    nameBodiesAndManoeuvresTrajectory.push_back("Venus");
    nameBodiesAndManoeuvresTrajectory.push_back("DSM3");
    nameBodiesAndManoeuvresTrajectory.push_back("Venus");
    nameBodiesAndManoeuvresTrajectory.push_back("DSM4");
    nameBodiesAndManoeuvresTrajectory.push_back("Mercury");

    std::vector< std::string > transferBodyTrajectory;
    transferBodyTrajectory.push_back("Earth");
    transferBodyTrajectory.push_back("Earth");
    transferBodyTrajectory.push_back("Venus");
    transferBodyTrajectory.push_back("Venus");
    transferBodyTrajectory.push_back("Mercury");


    std::vector< std::string > centralBody; centralBody.push_back( "Sun" );
    std::string bodyToPropagate = "spacecraft";


    // Create body map.
    simulation_setup::NamedBodyMap bodyMap = propagators::setupBodyMapFromEphemeridesForPatchedConicsTrajectory(centralBody[0],
            bodyToPropagate, transferBodyTrajectory);


    // Create acceleration map.
    std::vector< basic_astrodynamics::AccelerationMap > accelerationMap = propagators::setupAccelerationMapPatchedConicsTrajectory(
                transferBodyTrajectory.size(), centralBody[0], bodyToPropagate, bodyMap);


    // Create variable vector.
    std::vector< double > variableVector;

    // Add the time of flight and start epoch.
    variableVector.push_back( 1171.64503236 * physical_constants::JULIAN_DAY);
    variableVector.push_back( 399.999999715 * physical_constants::JULIAN_DAY);
    variableVector.push_back( 178.372255301 * physical_constants::JULIAN_DAY);
    variableVector.push_back( 299.223139512 * physical_constants::JULIAN_DAY);
    variableVector.push_back( 180.510754824 * physical_constants::JULIAN_DAY);
    variableVector.push_back( 1.0); // The capture time is irrelevant for the final leg.

    // Add the additional variables.
    // 1st leg.
    variableVector.push_back( 0.234594654679 );
    variableVector.push_back( 1408.99421278 );
    variableVector.push_back( 0.37992647165 * 2 * 3.14159265358979 );
    variableVector.push_back( std::acos(  2 * 0.498004040298 - 1. ) - 3.14159265358979 / 2 );
    // 2nd leg.
    variableVector.push_back( 0.0964769387134 );
    variableVector.push_back( 1.35077257078 );
    variableVector.push_back( 1.80629232251 * 6.378e6 );
    variableVector.push_back( 0.0 );
    // 3rd leg.
    variableVector.push_back( 0.829948744508);
    variableVector.push_back( 1.09554368115 );
    variableVector.push_back( 3.04129845698 * 6.052e6 );
    variableVector.push_back( 0.0 );
    // 4th leg.
    variableVector.push_back( 0.317174785637 );
    variableVector.push_back( 1.34317576594 );
    variableVector.push_back( 1.10000000891 * 6.052e6 );
    variableVector.push_back( 0.0 );



    // Create minimum pericenter radii vector
    std::vector< double > minimumPericenterRadii;
    minimumPericenterRadii.push_back( TUDAT_NAN ); minimumPericenterRadii.push_back( TUDAT_NAN ); minimumPericenterRadii.push_back( TUDAT_NAN );
    minimumPericenterRadii.push_back( TUDAT_NAN ); minimumPericenterRadii.push_back( TUDAT_NAN );

    // Create departure and capture variables.
    std::vector< double > semiMajorAxes;
    semiMajorAxes.push_back( std::numeric_limits< double >::infinity( ) ); semiMajorAxes.push_back( std::numeric_limits< double >::infinity( ) );
    std::vector< double > eccentricities;
    eccentricities.push_back( 0.0 ); eccentricities.push_back( 0.0 );


    // Define integrator settings.
    double initialTime = 0.0;
    double fixedStepSize = 1000.0;
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared < numerical_integrators::IntegratorSettings < > >(
                numerical_integrators::rungeKutta4, initialTime, fixedStepSize );



    // Compute difference between patched conics trajectory and full problem at departure and at arrival for each leg.
    // Compute difference between patched conics trajectory and full problem at departure and at arrival for each leg.
    for( int terminationType = 0; terminationType < 2; terminationType++ )
    {
        std::map< int, std::map< double, Eigen::Vector6d > > patchedConicsResultForEachLeg;
        std::map< int, std::map< double, Eigen::Vector6d > > fullProblemResultForEachLeg;
        std::map< int, std::map< double, Eigen::VectorXd > > dependentVariableResultForEachLeg;

        //    std::map< int, std::pair< Eigen::Vector6d, Eigen::Vector6d > > differenceStateArrivalAndDeparturePerLeg =
        tudat::propagators::fullPropagationPatchedConicsTrajectory(
                    bodyMap, accelerationMap, transferBodyTrajectory, centralBody[0], bodyToPropagate, legTypeVector, variableVector,
                minimumPericenterRadii, semiMajorAxes, eccentricities, integratorSettings, patchedConicsResultForEachLeg,
                fullProblemResultForEachLeg, dependentVariableResultForEachLeg,
                static_cast< bool >( terminationType ) );

        for( auto itr : patchedConicsResultForEachLeg )
        {
            for( auto innerItr : itr.second )
            {
                Eigen::Vector6d stateDifference = patchedConicsResultForEachLeg[ itr.first ][ innerItr.first ] -
                        fullProblemResultForEachLeg[ itr.first ][ innerItr.first ];
                for( int i = 0; i < 3; i++ )
                {
                    BOOST_CHECK_SMALL( std::fabs( stateDifference( i ) ), 1.0 );
                    BOOST_CHECK_SMALL( std::fabs( stateDifference( i + 3 ) ), 1.0E-6 );
                    BOOST_CHECK_SMALL( std::fabs( stateDifference( i ) ), 1.0 );
                    BOOST_CHECK_SMALL( std::fabs( stateDifference( i + 3 ) ), 1.0E-6 );
                }
            }
        }

        int numberOfLegs = patchedConicsResultForEachLeg.size( );
        auto legIterator = patchedConicsResultForEachLeg.begin( );
        auto nextLegIterator = patchedConicsResultForEachLeg.begin( );
        nextLegIterator++;

        for( int i = 0; i < numberOfLegs - 1; i++ )
        {
            double legTimeDifference =
                    nextLegIterator->second.begin( )->first - legIterator->second.rbegin( )->first;
            Eigen::Vector6d legStateDifference =
                    nextLegIterator->second.begin( )->second - legIterator->second.rbegin( )->second;
//            std::cout<<"Time diff. "<<legTimeDifference<<std::endl;
//            std::cout<<"State diff. "<<legStateDifference.transpose( )<<std::endl<<std::endl;

            if( terminationType == 0 )
            {
                BOOST_CHECK_SMALL( std::fabs( legTimeDifference ), 1.0E-6 );
                BOOST_CHECK_SMALL( std::fabs( legStateDifference.segment( 0, 3 ).norm( ) ), 1.0 );
            }
            else if( i % 2 == 0 )
            {
                BOOST_CHECK_SMALL( std::fabs( legTimeDifference ), 1.0E-6 );
                BOOST_CHECK_SMALL( std::fabs( legStateDifference.segment( 0, 3 ).norm( ) ), 1.0 );
            }
            else
            {
                BOOST_CHECK_EQUAL( std::fabs( legTimeDifference ) > 1.0E5, true );
            }

            legIterator++;
            nextLegIterator++;
        }
//        std::cout<<std::endl<<std::endl;
    }
}

}

}

}

