/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

//#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_MAIN

#include <limits>
#include <string>
#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/unitConversions.h"

#include <boost/test/unit_test.hpp>

#include <boost/lambda/lambda.hpp>

#include "tudat/astro/aerodynamics/exponentialAtmosphere.h"
#include "tudat/astro/basic_astro/sphericalStateConversions.h"
#include "tudat/astro/gravitation/centralGravityModel.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantDragCoefficient.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/empiricalAccelerationCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/gravitationalParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/radiationPressureCoefficient.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/ppnParameters.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/directTidalTimeLag.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantThrust.h"
#include "tudat/astro/relativity/einsteinInfeldHoffmannEquations.h"
#include "tudat/astro/orbit_determination/acceleration_partials/einsteinInfeldHoffmannPartials.h"
#include "tudat/astro/orbit_determination/acceleration_partials/numericalAccelerationPartial.h"
#include "tudat/astro/relativity/relativisticAccelerationCorrection.h"
#include "tudat/simulation/estimation_setup/createAccelerationPartials.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"
#include "tudat/simulation/estimation_setup/createEstimatableParameters.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/thrustSettings.h"
#include "tudat/simulation/environment_setup/createSystemModel.h"
//
//namespace tudat
//{
//namespace unit_tests
//{

using namespace tudat;
using namespace tudat::relativity;
using namespace tudat::gravitation;
using namespace tudat::aerodynamics;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::unit_conversions;
using namespace tudat::orbit_determination;
using namespace tudat::acceleration_partials;
using namespace tudat::spice_interface;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::electromagnetism;
using namespace tudat::basic_astrodynamics;
using namespace tudat::propulsion;

//BOOST_AUTO_TEST_SUITE( test_eih_partials )

//BOOST_AUTO_TEST_CASE( testEihPartials )
int main( )
{
    double testTime = 1.0E8;

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation end epoch.
    const double simulationStartEpoch = 0.0 * tudat::physical_constants::JULIAN_YEAR;
    const double simulationEndEpoch = 25.0 * tudat::physical_constants::JULIAN_YEAR;

    // Create body objects.
    std::vector<std::string> bodiesToCreate =
        { "Venus", "Sun" };
//        { "Earth", "Venus", "Mercury", "Sun" };//, "Moon", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune" };
    int numberOfBodies = bodiesToCreate.size( );
    BodyListSettings bodySettings =
        getDefaultBodySettings( bodiesToCreate );
    for( unsigned int i = 0; i < numberOfBodies; i++ )
    {
        bodySettings.at( bodiesToCreate.at( i ) )->gravityFieldSettings = centralGravityFromSpiceSettings( );
    }

    // Create Earth object
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    for( unsigned int i = 0; i < numberOfBodies; i++ )
    {
        bodies.at( bodiesToCreate.at( i ) )->setState(
            bodies.at( bodiesToCreate.at( i ) )->getStateInBaseFrameFromEphemeris( testTime ) );
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector<std::string> bodiesToPropagate = { "Venus" };//, "Venus" };
    std::vector<std::string> centralBodies = { "SSB" };

    int numberOfPropatedBodies = bodiesToPropagate.size( );

    // Define propagation settings.
    for ( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
    {
        std::map<std::string, std::vector<std::shared_ptr<AccelerationSettings> > > accelerationsOfCurrentBody;
        for ( unsigned int j = 0; j < numberOfBodies; j++ )
        {
            if ( bodiesToPropagate.at( i ) != bodiesToCreate.at( j ))
            {
                accelerationsOfCurrentBody[ bodiesToCreate.at( j ) ].push_back(
                    std::make_shared<AccelerationSettings>(
                        einstein_infeld_hoffmann_acceleration ) );
            }
        }
        accelerationMap[ bodiesToPropagate.at( i ) ] = accelerationsOfCurrentBody;
    }

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
        bodies, accelerationMap, bodiesToPropagate, centralBodies );

    std::shared_ptr< relativity::EinsteinInfeldHoffmannEquations > eihEquations = propagators::getEihEquationsFromAccelerationMap(
        accelerationModelMap );
    std::shared_ptr< acceleration_partials::EihEquationsPartials > eihPartials =
        std::make_shared< acceleration_partials::EihEquationsPartials >( eihEquations );
    eihPartials->update( testTime );

    double positionPerturbation = 1.0E6;
    double velocityPerturbation = 0.1;

    std::vector< std::vector< Eigen::Matrix< double, 1, 6 > > > upperturbedTotalPotentials_( numberOfBodies );
    std::vector< std::vector< Eigen::Matrix< double, 1, 6 > > > downperturbedTotalPotentials_( numberOfBodies );

    std::vector< std::vector< Eigen::Matrix< double, 1, 6 > > > upperturbedSingleContrubutionPotentials( numberOfBodies );
    std::vector< std::vector< Eigen::Matrix< double, 1, 6 > > > downperturbedSingleContrubutionPotentials( numberOfBodies );

    std::vector< std::vector< Eigen::Matrix< double, 3, 6 > > > upperturbedSinglePointmassAcceleration( numberOfBodies );
    std::vector< std::vector< Eigen::Matrix< double, 3, 6 > > > downperturbedSinglePointmassAcceleration( numberOfBodies );

    std::vector< std::vector< Eigen::Matrix< double, 3, 6 > > > upperturbedTotalPointmassAcceleration( numberOfBodies );
    std::vector< std::vector< Eigen::Matrix< double, 3, 6 > > > downperturbedTotalPointmassAcceleration( numberOfBodies );

    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 6 > > > > upperturbedScalarEihCorrections( 7 );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 6 > > > > downperturbedScalarEihCorrections( 7 );
    for( int k = 0; k < 7; k++ )
    {
        upperturbedScalarEihCorrections[ k ].resize( numberOfBodies );
        downperturbedScalarEihCorrections[ k ].resize( numberOfBodies );
    }

    std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > >  numericalTotalPotentialsWrtPosition_( numberOfBodies );
    std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > >  numericalSingleContributionPotentialsWrtPosition_( numberOfBodies );
    std::vector< std::vector< Eigen::Matrix< double, 3, 3 > > >  numericalSinglePointMassAccelerationWrtPosition_( numberOfBodies );
    std::vector< std::vector< Eigen::Matrix< double, 3, 3 > > >  numericalTotalPointMassAccelerationWrtPosition_( numberOfBodies );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > > > numericalScalarEihCorrectionsWrtExertingPosition_( 7 );

    std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > >  analyticalTotalPotentialsWrtPosition_ = eihPartials->getCurrentTotalPotentialWrtPosition( );
    std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > >  analyticalSingleContributionPotentialsWrtPosition_ = eihPartials->getCurrentLocalPotentialWrtPosition( );
    std::vector< std::vector< Eigen::Matrix< double, 3, 3 > > >  analyticalSinglePointMassAccelerationWrtPosition_ = eihPartials->getCurrentSinglePointMassAccelerationWrtPosition( );
    std::vector< std::vector< Eigen::Matrix< double, 3, 3 > > >  analyticalTotalPointMassAccelerationWrtPosition_ = eihPartials->getCurrentTotalPointMassAccelerationWrtPosition( );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > > > analyticalScalarEihCorrectionsWrtExertingPosition_( 7 );

    for( int k = 0; k < 7; k++ )
    {
        numericalScalarEihCorrectionsWrtExertingPosition_[ k ].resize( numberOfBodies );
        analyticalScalarEihCorrectionsWrtExertingPosition_[ k ].resize( numberOfBodies );
    }
    
    for( unsigned int i = 0; i < numberOfBodies; i++ )
    {
        upperturbedTotalPotentials_[ i ].resize( numberOfBodies );
        downperturbedTotalPotentials_[ i ].resize( numberOfBodies );

        upperturbedSingleContrubutionPotentials[ i ].resize( numberOfBodies );
        downperturbedSingleContrubutionPotentials[ i ].resize( numberOfBodies );

        upperturbedSinglePointmassAcceleration[ i ].resize( numberOfBodies );
        downperturbedSinglePointmassAcceleration[ i ].resize( numberOfBodies );

        upperturbedTotalPointmassAcceleration[ i ].resize( numberOfBodies );
        downperturbedTotalPointmassAcceleration[ i ].resize( numberOfBodies );

        numericalTotalPotentialsWrtPosition_[ i ].resize( numberOfBodies );
        numericalSingleContributionPotentialsWrtPosition_[ i ].resize( numberOfBodies );
        numericalSinglePointMassAccelerationWrtPosition_[ i ].resize( numberOfBodies );
        numericalTotalPointMassAccelerationWrtPosition_[ i ].resize( numberOfBodies );

        for( int k = 0; k < 7; k++ )
        {
            upperturbedScalarEihCorrections[ k ][ i ].resize( numberOfBodies );
            downperturbedScalarEihCorrections[ k ][ i ].resize( numberOfBodies );
            numericalScalarEihCorrectionsWrtExertingPosition_[ k ][ i ].resize( numberOfBodies );
            analyticalScalarEihCorrectionsWrtExertingPosition_[ k ][ i ].resize( numberOfBodies );
        }

    }
    for( unsigned int i = 0; i < numberOfBodies; i++ )
    {
        Eigen::Vector6d nominalState = bodies.at( bodiesToCreate.at( i ) )->getState( );
        for( unsigned int index = 0; index < 6; index++ )
        {
            eihEquations->update( TUDAT_NAN );
            eihEquations->update( testTime );
            for( int k = 0; k < 7; k++ )
            {
                for( unsigned int j = 0; j < numberOfPropatedBodies; j++ )
                {
                    analyticalScalarEihCorrectionsWrtExertingPosition_[ k ][ j ][ i ].setZero( );
                    eihPartials->addSingleScalarTermWrtPositionPartial(
                        analyticalScalarEihCorrectionsWrtExertingPosition_[ k ][ j ][ i ], j, i, true, k );
                }
            }
            Eigen::Vector6d upperturbedState = nominalState;
            upperturbedState( index ) += ( ( index < 3 ) ? positionPerturbation : velocityPerturbation );
            bodies.at( bodiesToCreate.at( i ) )->setState( upperturbedState );
            eihEquations->update( TUDAT_NAN );
            eihEquations->update( testTime );

            for( unsigned int j = 0; j < numberOfPropatedBodies; j++ )
            {
                upperturbedTotalPotentials_[ i ][ j ]( index ) = eihEquations->getLocalPotential( j );
                upperturbedSingleContrubutionPotentials[ i ][ j ]( index ) = eihEquations->getSingleSourceLocalPotential( j, i );
                upperturbedSinglePointmassAcceleration[ i ][ j ].block( 0, index, 3, 1 )= eihEquations->getSinglePointMassAccelerations( j, i );
                upperturbedTotalPointmassAcceleration[ i ][ j ].block( 0, index, 3, 1 ) = eihEquations->getTotalPointMassAcceleration( j );
                for( int k = 0; k < 7; k++ )
                {
                    upperturbedScalarEihCorrections[ k ][ i ][ j ]( index ) = eihEquations->getScalarEihCorrections( ).at( k ).at( j ).at( i );
                }
            }

            Eigen::Vector6d downperturbedState = nominalState;
            downperturbedState( index ) -= ( ( index < 3 ) ? positionPerturbation : velocityPerturbation );
            bodies.at( bodiesToCreate.at( i ) )->setState( downperturbedState );
            eihEquations->update( TUDAT_NAN );
            eihEquations->update( testTime );
            for( unsigned int j = 0; j < numberOfPropatedBodies; j++ )
            {
                downperturbedTotalPotentials_[ i ][ j ]( index ) = eihEquations->getLocalPotential( j );
                downperturbedSingleContrubutionPotentials[ i ][ j ]( index ) = eihEquations->getSingleSourceLocalPotential( j, i );
                downperturbedSinglePointmassAcceleration[ i ][ j ].block( 0, index, 3, 1 )= eihEquations->getSinglePointMassAccelerations( j, i );
                downperturbedTotalPointmassAcceleration[ i ][ j ].block( 0, index, 3, 1 ) = eihEquations->getTotalPointMassAcceleration( j );
                for( int k = 0; k < 7; k++ )
                {
                    downperturbedScalarEihCorrections[ k ][ i ][ j ]( index ) = eihEquations->getScalarEihCorrections( ).at( k ).at( j ).at( i );
                }
            }

            bodies.at( bodiesToCreate.at( i ) )->setState( nominalState );

        }

        for( int j = 0; j < numberOfPropatedBodies; j++ )
        {
            {
                numericalTotalPotentialsWrtPosition_[ i ][ j ] =
                    ( upperturbedTotalPotentials_[ i ][ j ].block( 0, 0, 1, 3 ) - downperturbedTotalPotentials_[ i ][ j ].block( 0, 0, 1, 3 ) ) /
                    ( 2.0 * positionPerturbation );
                numericalSingleContributionPotentialsWrtPosition_[ i ][ j ] =
                    ( upperturbedSingleContrubutionPotentials[ i ][ j ].block( 0, 0, 1, 3 ) - downperturbedSingleContrubutionPotentials[ i ][ j ].block( 0, 0, 1, 3 ) ) /
                    ( 2.0 * positionPerturbation );
                numericalSinglePointMassAccelerationWrtPosition_[ i ][ j ] =
                    ( upperturbedSinglePointmassAcceleration[ i ][ j ].block( 0, 0, 3, 3 ) - downperturbedSinglePointmassAcceleration[ i ][ j ].block( 0, 0, 3, 3 ) ) /
                    ( 2.0 * positionPerturbation );
                numericalTotalPointMassAccelerationWrtPosition_[ i ][ j ] =
                    ( upperturbedTotalPointmassAcceleration[ i ][ j ].block( 0, 0, 3, 3 ) - downperturbedTotalPointmassAcceleration[ i ][ j ].block( 0, 0, 3, 3 ) ) /
                    ( 2.0 * positionPerturbation );
                if( i != j )
                {
                    for( int k = 0; k < 7; k++ )
                    {
                        numericalScalarEihCorrectionsWrtExertingPosition_[ k ][ i ][ j ] =
                            ( upperturbedScalarEihCorrections[ k ][ i ][ j ].block( 0, 0, 1, 3 ) - downperturbedScalarEihCorrections[ k ][ i ][ j ].block( 0, 0, 1, 3 ) ) /
                            ( 2.0 * positionPerturbation );
                        std::cout<<k<<" "<<i<<" "<<j<<std::endl<<
                        upperturbedScalarEihCorrections[ k ][ i ][ j ]<<std::endl<<
                        downperturbedScalarEihCorrections[ k ][ i ][ j ]<<std::endl<<
                        numericalScalarEihCorrectionsWrtExertingPosition_[ k ][ i ][ j ]<<std::endl<<
                        analyticalScalarEihCorrectionsWrtExertingPosition_[ k ][ j ][ i ]<<std::endl<<std::endl;
                    }
                }

//                std::cout<<( analyticalTotalPotentialsWrtPosition_[ j ][ i ].block( 0, 0, 1, 3 ) -
//                    numericalTotalPotentialsWrtPosition_[ i ][ j ].block( 0, 0, 1, 3 ) ).
//                    cwiseQuotient( numericalTotalPotentialsWrtPosition_[ i ][ j ].block( 0, 0, 1, 3 ) )<<std::endl<<std::endl;
//
//                std::cout<<( analyticalTotalPointMassAccelerationWrtPosition_[ j ][ i ].block( 0, 0, 3, 3 ) -
//                             numericalTotalPointMassAccelerationWrtPosition_[ i ][ j ].block( 0, 0, 3, 3 ) ).
//                    cwiseQuotient( numericalTotalPointMassAccelerationWrtPosition_[ i ][ j ].block( 0, 0, 3, 3 ) )<<std::endl<<std::endl;
//
//                if( i != j )
//                {
//                    std::cout<<( analyticalSingleContributionPotentialsWrtPosition_[ j ][ i ].block( 0, 0, 1, 3 ) -
//                        numericalSingleContributionPotentialsWrtPosition_[ i ][ j ].block( 0, 0, 1, 3 ) ).
//                        cwiseQuotient( numericalSingleContributionPotentialsWrtPosition_[ i ][ j ].block( 0, 0, 1, 3 ) )<<std::endl<<std::endl;
//
//                    std::cout<<( analyticalSinglePointMassAccelerationWrtPosition_[ j ][ i ].block( 0, 0, 3, 3 ) -
//                        numericalSinglePointMassAccelerationWrtPosition_[ i ][ j ].block( 0, 0, 3, 3 ) ).
//                        cwiseQuotient( numericalSinglePointMassAccelerationWrtPosition_[ i ][ j ].block( 0, 0, 3, 3 ) )<<std::endl<<std::endl;
//                }

            }

        }
    }


}

//BOOST_AUTO_TEST_SUITE_END( )
//
//} // namespace unit_tests
//
//} // namespace tudat
