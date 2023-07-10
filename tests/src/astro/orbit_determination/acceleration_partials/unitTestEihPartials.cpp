/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

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

namespace tudat
{
namespace unit_tests
{

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

BOOST_AUTO_TEST_SUITE( test_eih_partials )

BOOST_AUTO_TEST_CASE( testEihPartials )
//int main( )
{
    double testTime = 1.0E8;

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation end epoch.
    const double simulationStartEpoch = 0.0 * tudat::physical_constants::JULIAN_YEAR;
    const double simulationEndEpoch = 25.0 * tudat::physical_constants::JULIAN_YEAR;

    // Create body objects.
    std::vector<std::string> bodiesToCreate =
        { "Jupiter", "Sun" , "Venus" };
//        { "Earth", "Venus", "Mercury", "Sun" };//, "Moon", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune" };
    unsigned int numberOfBodies = bodiesToCreate.size( );
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
    std::vector<std::string> bodiesToPropagate = { "Jupiter", "Sun", "Venus" };//, "Venus" };
    std::vector<std::string> centralBodies = { "SSB", "SSB", "SSB" };

    unsigned int numberOfPropatedBodies = bodiesToPropagate.size( );

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

    double positionPerturbation = 1.0E8;
    double velocityPerturbation = 0.1;

    std::vector< std::vector< Eigen::Matrix< double, 1, 6 > > > upperturbedTotalPotentials_( numberOfBodies );
    std::vector< std::vector< Eigen::Matrix< double, 1, 6 > > > downperturbedTotalPotentials_( numberOfBodies );

    std::vector< std::vector< Eigen::Matrix< double, 1, 6 > > > upperturbedSingleContrubutionPotentials( numberOfBodies );
    std::vector< std::vector< Eigen::Matrix< double, 1, 6 > > > downperturbedSingleContrubutionPotentials( numberOfBodies );

    std::vector< std::vector< Eigen::Matrix< double, 3, 6 > > > upperturbedSinglePointmassAcceleration( numberOfBodies );
    std::vector< std::vector< Eigen::Matrix< double, 3, 6 > > > downperturbedSinglePointmassAcceleration( numberOfBodies );

    std::vector< std::vector< Eigen::Matrix< double, 3, 6 > > > upperturbedTotalPointmassAcceleration( numberOfBodies );
    std::vector< std::vector< Eigen::Matrix< double, 3, 6 > > > downperturbedTotalPointmassAcceleration( numberOfBodies );

    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 6 > > > > upperturbedExertingScalarEihCorrections( 7 );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 6 > > > > downperturbedExertingScalarEihCorrections( 7 );

    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 6 > > > > upperturbedUndergoingScalarEihCorrections( 7 );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 6 > > > > downperturbedUndergoingScalarEihCorrections( 7 );

    std::vector< std::vector< std::vector< Eigen::Matrix< double, 3, 6 > > > > upperturbedExertingVectorEihCorrections( 3 );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 3, 6 > > > > downperturbedExertingVectorEihCorrections( 3 );

    std::vector< std::vector< std::vector< Eigen::Matrix< double, 3, 6 > > > > upperturbedUndergoingVectorEihCorrections( 3 );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 3, 6 > > > > downperturbedUndergoingVectorEihCorrections( 3 );

    for( int k = 0; k < 7; k++ )
    {
        upperturbedExertingScalarEihCorrections[ k ].resize( numberOfBodies );
        downperturbedExertingScalarEihCorrections[ k ].resize( numberOfBodies );

        upperturbedUndergoingScalarEihCorrections[ k ].resize( numberOfBodies );
        downperturbedUndergoingScalarEihCorrections[ k ].resize( numberOfBodies );

        if( k < 3 )
        {
            upperturbedExertingVectorEihCorrections[ k ].resize( numberOfBodies );
            downperturbedExertingVectorEihCorrections[ k ].resize( numberOfBodies );

            upperturbedUndergoingVectorEihCorrections[ k ].resize( numberOfBodies );
            downperturbedUndergoingVectorEihCorrections[ k ].resize( numberOfBodies );
        }
    }

    std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > >  numericalTotalPotentialsWrtPosition_( numberOfBodies );
    std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > >  numericalSingleContributionPotentialsWrtPosition_( numberOfBodies );
    std::vector< std::vector< Eigen::Matrix< double, 3, 3 > > >  numericalSinglePointMassAccelerationWrtPosition_( numberOfBodies );
    std::vector< std::vector< Eigen::Matrix< double, 3, 3 > > >  numericalTotalPointMassAccelerationWrtPosition_( numberOfBodies );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > > > numericalScalarEihCorrectionsWrtExertingPosition_( 7 );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > > > numericalScalarEihCorrectionsWrtExertingVelocity_( 7 );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > > > numericalScalarEihCorrectionsWrtUndergoingPosition_( 7 );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > > > numericalVectorEihCorrectionsWrtExertingPosition_( 3 );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > > > numericalVectorEihCorrectionsWrtExertingVelocity_( 3 );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > > > numericalVectorEihCorrectionsWrtUndergoingPosition_( 3 );

    std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > >  analyticalTotalPotentialsWrtPosition_ = eihPartials->getCurrentTotalPotentialWrtPosition( );
    std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > >  analyticalSingleContributionPotentialsWrtPosition_ = eihPartials->getCurrentLocalPotentialWrtPosition( );
    std::vector< std::vector< Eigen::Matrix< double, 3, 3 > > >  analyticalSinglePointMassAccelerationWrtPosition_ = eihPartials->getCurrentSinglePointMassAccelerationWrtPosition( );
    std::vector< std::vector< Eigen::Matrix< double, 3, 3 > > >  analyticalTotalPointMassAccelerationWrtPosition_ = eihPartials->getCurrentTotalPointMassAccelerationWrtPosition( );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > > > analyticalScalarEihCorrectionsWrtExertingPosition_( 7 );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > > > analyticalScalarEihCorrectionsWrtExertingVelocity_( 7 );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > > > analyticalScalarEihCorrectionsWrtUndergoingPosition_( 7 );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 3, 3 > > > > analyticalVectorEihCorrectionsWrtExertingPosition_( 3 );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 3, 3 > > > > analyticalVectorEihCorrectionsWrtExertingVelocity_( 3 );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 3, 3 > > > > analyticalVectorEihCorrectionsWrtUndergoingPosition_( 3 );

    for( int k = 0; k < 7; k++ )
    {
        numericalScalarEihCorrectionsWrtExertingPosition_[ k ].resize( numberOfBodies );
        analyticalScalarEihCorrectionsWrtExertingPosition_[ k ].resize( numberOfBodies );

        numericalScalarEihCorrectionsWrtExertingVelocity_[ k ].resize( numberOfBodies );
        analyticalScalarEihCorrectionsWrtExertingVelocity_[ k ].resize( numberOfBodies );

        numericalScalarEihCorrectionsWrtUndergoingPosition_[ k ].resize( numberOfBodies );
        analyticalScalarEihCorrectionsWrtUndergoingPosition_[ k ].resize( numberOfBodies );

        if( k < 3 )
        {
            numericalVectorEihCorrectionsWrtExertingPosition_[ k ].resize( numberOfBodies );
            analyticalVectorEihCorrectionsWrtExertingPosition_[ k ].resize( numberOfBodies );

            numericalVectorEihCorrectionsWrtExertingVelocity_[ k ].resize( numberOfBodies );
            analyticalVectorEihCorrectionsWrtExertingVelocity_[ k ].resize( numberOfBodies );

            numericalVectorEihCorrectionsWrtUndergoingPosition_[ k ].resize( numberOfBodies );
            analyticalVectorEihCorrectionsWrtUndergoingPosition_[ k ].resize( numberOfBodies );
        }
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
            upperturbedExertingScalarEihCorrections[ k ][ i ].resize( numberOfBodies );
            downperturbedExertingScalarEihCorrections[ k ][ i ].resize( numberOfBodies );

            upperturbedUndergoingScalarEihCorrections[ k ][ i ].resize( numberOfBodies );
            downperturbedUndergoingScalarEihCorrections[ k ][ i ].resize( numberOfBodies );

            numericalScalarEihCorrectionsWrtExertingPosition_[ k ][ i ].resize( numberOfBodies );
            analyticalScalarEihCorrectionsWrtExertingPosition_[ k ][ i ].resize( numberOfBodies );

            numericalScalarEihCorrectionsWrtExertingVelocity_[ k ][ i ].resize( numberOfBodies );
            analyticalScalarEihCorrectionsWrtExertingVelocity_[ k ][ i ].resize( numberOfBodies );

            numericalScalarEihCorrectionsWrtUndergoingPosition_[ k ][ i ].resize( numberOfBodies );
            analyticalScalarEihCorrectionsWrtUndergoingPosition_[ k ][ i ].resize( numberOfBodies );

            if( k < 3 )
            {
                upperturbedExertingVectorEihCorrections[ k ][ i ].resize( numberOfBodies );
                downperturbedExertingVectorEihCorrections[ k ][ i ].resize( numberOfBodies );

                upperturbedUndergoingVectorEihCorrections[ k ][ i ].resize( numberOfBodies );
                downperturbedUndergoingVectorEihCorrections[ k ][ i ].resize( numberOfBodies );

                numericalVectorEihCorrectionsWrtExertingPosition_[ k ][ i ].resize( numberOfBodies );
                analyticalVectorEihCorrectionsWrtExertingPosition_[ k ][ i ].resize( numberOfBodies );

                numericalVectorEihCorrectionsWrtExertingVelocity_[ k ][ i ].resize( numberOfBodies );
                analyticalVectorEihCorrectionsWrtExertingVelocity_[ k ][ i ].resize( numberOfBodies );

                numericalVectorEihCorrectionsWrtUndergoingPosition_[ k ][ i ].resize( numberOfBodies );
                analyticalVectorEihCorrectionsWrtUndergoingPosition_[ k ][ i ].resize( numberOfBodies );
            }
        }

    }
    for( unsigned int i = 0; i < numberOfBodies; i++ )
    {
        // Get current body nominal state
        Eigen::Vector6d nominalState = bodies.at( bodiesToCreate.at( i ))->getState( );
        for ( unsigned int index = 0; index < 6; index++ )
        {
            // Reset to nominal state
            bodies.at( bodiesToCreate.at( i ))->setState( nominalState );
            eihEquations->update( TUDAT_NAN );
            eihEquations->update( testTime );

            // Iterate over each scalar term
            for ( int k = 0; k < 7; k++ )
            {
                for ( unsigned int j = 0; j < numberOfPropatedBodies; j++ )
                {
                    // Set d(fs)_{k,ji}/dx_{i} to 0
                    analyticalScalarEihCorrectionsWrtExertingPosition_[ k ][ j ][ i ].setZero( );

                    analyticalScalarEihCorrectionsWrtExertingVelocity_[ k ][ j ][ i ].setZero( );

                    // Set d(fs)_{k,ji}/dx_{j} to 0
                    analyticalScalarEihCorrectionsWrtUndergoingPosition_[ k ][ j ][ i ].setZero( );

                    // Compute analytical d(fs)_{k,ji}/dx_{i}
                    eihPartials->addSingleScalarTermWrtPositionPartial(
                        analyticalScalarEihCorrectionsWrtExertingPosition_[ k ][ j ][ i ], j, i, true, k );
                    eihPartials->addSingleScalarTermWrtVelocityPartial(
                        analyticalScalarEihCorrectionsWrtExertingVelocity_[ k ][ j ][ i ], j, i, true, k );

                    // Compute analytical d(fs)_{k,ji}/dx_{j}
                    eihPartials->addSingleScalarTermWrtPositionPartial(
                        analyticalScalarEihCorrectionsWrtUndergoingPosition_[ k ][ j ][ i ], j, i, false, k );

                    if ( k < 3 )
                    {
                        // Set d(fv)_{k,ji}/dx_{i} to 0
                        analyticalVectorEihCorrectionsWrtExertingPosition_[ k ][ j ][ i ].setZero( );

                        analyticalVectorEihCorrectionsWrtExertingVelocity_[ k ][ j ][ i ].setZero( );

                        // Set d(fv)_{k,ji}/dx_{j} to 0
                        analyticalVectorEihCorrectionsWrtUndergoingPosition_[ k ][ j ][ i ].setZero( );

                        // Compute analytical d(fv)_{k,ji}/dx_{i}
                        eihPartials->addSingleVectorTermWrtPositionPartial(
                            analyticalVectorEihCorrectionsWrtExertingPosition_[ k ][ j ][ i ], j, i, true, k );

                        eihPartials->addSingleVectorTermWrtVelocityPartial(
                            analyticalVectorEihCorrectionsWrtExertingVelocity_[ k ][ j ][ i ], j, i, true, k );

                        // Compute analytical d(fv)_{k,ji}/dx_{j}
                        eihPartials->addSingleVectorTermWrtPositionPartial(
                            analyticalVectorEihCorrectionsWrtUndergoingPosition_[ k ][ j ][ i ], j, i, false, k );
                    }
                }
            }

            // Perturb exerting body (x_{i}) up
            Eigen::Vector6d upperturbedState = nominalState;
            upperturbedState( index ) += (( index < 3 ) ? positionPerturbation : velocityPerturbation );
            bodies.at( bodiesToCreate.at( i ))->setState( upperturbedState );
            eihEquations->update( TUDAT_NAN );
            eihEquations->update( testTime );

            for ( unsigned int j = 0; j < numberOfPropatedBodies; j++ )
            {
                // Compute upperturbed potentials and accelerations
                upperturbedTotalPotentials_[ i ][ j ]( index ) = eihEquations->getLocalPotential( j );
                upperturbedSingleContrubutionPotentials[ i ][ j ]( index ) =
                    eihEquations->getSingleSourceLocalPotential( j, i );
                upperturbedSinglePointmassAcceleration[ i ][ j ].block( 0, index, 3, 1 ) =
                    eihEquations->getSinglePointMassAccelerations( j, i );
                upperturbedTotalPointmassAcceleration[ i ][ j ].block( 0, index, 3, 1 ) =
                    eihEquations->getTotalPointMassAcceleration( j );

                // Compute upperturbed scalar and vector terms
                for ( int k = 0; k < 7; k++ )
                {
                    upperturbedExertingScalarEihCorrections[ k ][ i ][ j ]( index ) =
                        eihEquations->getScalarEihCorrections( ).at( k ).at( j ).at( i );
                    if ( j < numberOfPropatedBodies )
                    {
                        upperturbedUndergoingScalarEihCorrections[ k ][ i ][ j ]( index ) =
                            eihEquations->getScalarEihCorrections( ).at( k ).at( j ).at( i );
                    }

                    if ( k < 3 )
                    {
                        upperturbedExertingVectorEihCorrections[ k ][ i ][ j ].block( 0, index, 3, 1 ) =
                            eihEquations->getVectorEihCorrections( ).at( k ).at( j ).at( i );
                        if ( j < numberOfPropatedBodies )
                        {
                            upperturbedUndergoingVectorEihCorrections[ k ][ i ][ j ].block( 0, index, 3, 1 ) =
                                eihEquations->getVectorEihCorrections( ).at( k ).at( i ).at( j );
                        }
                    }
                }

            }

            Eigen::Vector6d downperturbedState = nominalState;
            downperturbedState( index ) -= (( index < 3 ) ? positionPerturbation : velocityPerturbation );
            bodies.at( bodiesToCreate.at( i ))->setState( downperturbedState );
            eihEquations->update( TUDAT_NAN );
            eihEquations->update( testTime );
            for ( unsigned int j = 0; j < numberOfPropatedBodies; j++ )
            {
                downperturbedTotalPotentials_[ i ][ j ]( index ) = eihEquations->getLocalPotential( j );
                downperturbedSingleContrubutionPotentials[ i ][ j ]( index ) =
                    eihEquations->getSingleSourceLocalPotential( j, i );
                downperturbedSinglePointmassAcceleration[ i ][ j ].block( 0, index, 3, 1 ) =
                    eihEquations->getSinglePointMassAccelerations( j, i );
                downperturbedTotalPointmassAcceleration[ i ][ j ].block( 0, index, 3, 1 ) =
                    eihEquations->getTotalPointMassAcceleration( j );
                for ( int k = 0; k < 7; k++ )
                {
                    downperturbedExertingScalarEihCorrections[ k ][ i ][ j ]( index ) =
                        eihEquations->getScalarEihCorrections( ).at( k ).at( j ).at( i );
                    if ( j < numberOfPropatedBodies )
                    {
                        downperturbedUndergoingScalarEihCorrections[ k ][ i ][ j ]( index ) =
                            eihEquations->getScalarEihCorrections( ).at( k ).at( j ).at( i );
                    }

                    if ( k < 3 )
                    {
                        downperturbedExertingVectorEihCorrections[ k ][ i ][ j ].block( 0, index, 3, 1 ) =
                            eihEquations->getVectorEihCorrections( ).at( k ).at( j ).at( i );
                        if ( j < numberOfPropatedBodies )
                        {
                            downperturbedUndergoingVectorEihCorrections[ k ][ i ][ j ].block( 0, index, 3, 1 ) =
                                eihEquations->getVectorEihCorrections( ).at( k ).at( i ).at( j );
                        }
                    }
                }
            }
        }
    }

    for( unsigned int i = 0; i < numberOfBodies; i++ )
    {
        for( unsigned int j = 0; j < numberOfPropatedBodies; j++ )
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
                        ( upperturbedExertingScalarEihCorrections[ k ][ i ][ j ].block( 0, 0, 1, 3 ) - downperturbedExertingScalarEihCorrections[ k ][ i ][ j ].block( 0, 0, 1, 3 ) ) /
                        ( 2.0 * positionPerturbation );

                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        numericalScalarEihCorrectionsWrtExertingPosition_[ k ][ i ][ j ].block( 0, 0, 1, 3 ),
                        analyticalScalarEihCorrectionsWrtExertingPosition_[ k ][ j ][ i ].block( 0, 0, 1, 3 ),
                        1.0E-4);

                    numericalScalarEihCorrectionsWrtExertingVelocity_[ k ][ i ][ j ] =
                        ( upperturbedExertingScalarEihCorrections[ k ][ i ][ j ].block( 0, 3, 1, 3 ) - downperturbedExertingScalarEihCorrections[ k ][ i ][ j ].block( 0, 3, 1, 3 ) ) /
                        ( 2.0 * velocityPerturbation );

                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        numericalScalarEihCorrectionsWrtExertingVelocity_[ k ][ i ][ j ].block( 0, 0, 1, 3 ),
                        analyticalScalarEihCorrectionsWrtExertingVelocity_[ k ][ j ][ i ].block( 0, 0, 1, 3 ),
                        1.0E-4);

//                    std::cout<<"k "<<k<<" (i,j) "<<i<<" "<<j<<std::endl
//                             <<numericalScalarEihCorrectionsWrtExertingVelocity_[ k ][ i ][ j ]<<std::endl
//                             <<analyticalScalarEihCorrectionsWrtExertingVelocity_[ k ][ j ][ i ]<<std::endl<<std::endl;

                }
            }

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                numericalTotalPotentialsWrtPosition_[ i ][ j ].block( 0, 0, 1, 3 ),
                analyticalTotalPotentialsWrtPosition_[ j ][ i ].block( 0, 0, 1, 3 ),
                1.0E-4);

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                numericalTotalPointMassAccelerationWrtPosition_[ i ][ j ].block( 0, 0, 3, 3 ),
                analyticalTotalPointMassAccelerationWrtPosition_[ j ][ i ].block( 0, 0, 3, 3 ),
                1.0E-4);

            if( i != j )
            {

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    numericalSingleContributionPotentialsWrtPosition_[ i ][ j ].block( 0, 0, 1, 3 ),
                    analyticalSingleContributionPotentialsWrtPosition_[ j ][ i ].block( 0, 0, 1, 3 ),
                    1.0E-4);

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    numericalSinglePointMassAccelerationWrtPosition_[ i ][ j ].block( 0, 0, 3, 3 ),
                    analyticalSinglePointMassAccelerationWrtPosition_[ j ][ i ].block( 0, 0, 3, 3 ),
                    1.0E-4);
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
