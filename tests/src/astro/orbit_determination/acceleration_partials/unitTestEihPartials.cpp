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

    std::vector< std::vector< Eigen::Matrix< double, 1, 6 > > > upperturbedTotalPotentials_ = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 1, 6 > >(
        numberOfBodies, numberOfBodies, Eigen::Matrix< double, 1, 6 >::Zero( ) );
    std::vector< std::vector< Eigen::Matrix< double, 1, 6 > > > downperturbedTotalPotentials_ = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 1, 6 > >(
        numberOfBodies, numberOfBodies, Eigen::Matrix< double, 1, 6 >::Zero( ) );

    std::vector< std::vector< Eigen::Matrix< double, 1, 6 > > > upperturbedSingleContrubutionPotentials = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 1, 6 > >(
        numberOfBodies, numberOfBodies, Eigen::Matrix< double, 1, 6 >::Zero( ) );
    std::vector< std::vector< Eigen::Matrix< double, 1, 6 > > > downperturbedSingleContrubutionPotentials = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 1, 6 > >(
        numberOfBodies, numberOfBodies, Eigen::Matrix< double, 1, 6 >::Zero( ) );

    std::vector< std::vector< Eigen::Matrix< double, 3, 6 > > > upperturbedSinglePointmassAcceleration = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 3, 6 > >(
        numberOfBodies, numberOfBodies, Eigen::Matrix< double, 3, 6 >::Zero( ) );
    std::vector< std::vector< Eigen::Matrix< double, 3, 6 > > > downperturbedSinglePointmassAcceleration = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 3, 6 > >(
        numberOfBodies, numberOfBodies, Eigen::Matrix< double, 3, 6 >::Zero( ) );

    std::vector< std::vector< Eigen::Matrix< double, 3, 6 > > > upperturbedTotalPointmassAcceleration = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 3, 6 > >(
        numberOfBodies, numberOfBodies, Eigen::Matrix< double, 3, 6 >::Zero( ) );
    std::vector< std::vector< Eigen::Matrix< double, 3, 6 > > > downperturbedTotalPointmassAcceleration = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 3, 6 > >(
        numberOfBodies, numberOfBodies, Eigen::Matrix< double, 3, 6 >::Zero( ) );

    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 6 > > > > upperturbedExertingScalarEihCorrections = utilities::getThreeDimensionalVector< Eigen::Matrix< double, 1, 6 > >(
        7, numberOfBodies, numberOfBodies, Eigen::Matrix< double, 1, 6 >::Zero( ) );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 6 > > > > downperturbedExertingScalarEihCorrections = utilities::getThreeDimensionalVector< Eigen::Matrix< double, 1, 6 > >(
        7, numberOfBodies, numberOfBodies, Eigen::Matrix< double, 1, 6 >::Zero( ) );

    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 6 > > > > upperturbedUndergoingScalarEihCorrections = utilities::getThreeDimensionalVector< Eigen::Matrix< double, 1, 6 > >(
        7, numberOfBodies, numberOfBodies, Eigen::Matrix< double, 1, 6 >::Zero( ) );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 6 > > > > downperturbedUndergoingScalarEihCorrections = utilities::getThreeDimensionalVector< Eigen::Matrix< double, 1, 6 > >(
        7, numberOfBodies, numberOfBodies, Eigen::Matrix< double, 1, 6 >::Zero( ) );

    std::vector< std::vector< std::vector< Eigen::Matrix< double, 3, 6 > > > > upperturbedExertingVectorEihCorrections = utilities::getThreeDimensionalVector< Eigen::Matrix< double, 3, 6 > >(
        3, numberOfBodies, numberOfBodies, Eigen::Matrix< double, 3, 6 >::Zero( ) );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 3, 6 > > > > downperturbedExertingVectorEihCorrections = utilities::getThreeDimensionalVector< Eigen::Matrix< double, 3, 6 > >(
        3, numberOfBodies, numberOfBodies, Eigen::Matrix< double, 3, 6 >::Zero( ) );

    std::vector< std::vector< std::vector< Eigen::Matrix< double, 3, 6 > > > > upperturbedUndergoingVectorEihCorrections = utilities::getThreeDimensionalVector< Eigen::Matrix< double, 3, 6 > >(
        3, numberOfBodies, numberOfBodies, Eigen::Matrix< double, 3, 6 >::Zero( ) );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 3, 6 > > > > downperturbedUndergoingVectorEihCorrections = utilities::getThreeDimensionalVector< Eigen::Matrix< double, 3, 6 > >(
        3, numberOfBodies, numberOfBodies, Eigen::Matrix< double, 3, 6 >::Zero( ) );


    std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > >  numericalTotalPotentialsWrtPosition = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 1, 3 > >(
        numberOfBodies, numberOfBodies, Eigen::Matrix< double, 1, 3 >::Zero( ) );
    std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > >  numericalSingleContributionPotentialsWrtPosition = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 1, 3 > >(
        numberOfBodies, numberOfBodies, Eigen::Matrix< double, 1, 3 >::Zero( ) );
    std::vector< std::vector< Eigen::Matrix< double, 3, 3 > > >  numericalSinglePointMassAccelerationWrtPosition = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 3, 3 > >(
        numberOfBodies, numberOfBodies, Eigen::Matrix< double, 3, 3 >::Zero( ) );
    std::vector< std::vector< Eigen::Matrix< double, 3, 3 > > >  numericalTotalPointMassAccelerationWrtPosition = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 3, 3 > >(
        numberOfBodies, numberOfBodies, Eigen::Matrix< double, 3, 3 >::Zero( ) );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > > > numericalScalarEihCorrectionsWrtExertingPosition = utilities::getThreeDimensionalVector< Eigen::Matrix< double, 1, 3 > >(
        7, numberOfBodies, numberOfBodies, Eigen::Matrix< double, 1, 3 >::Zero( ) );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > > > numericalScalarEihCorrectionsWrtExertingVelocity = utilities::getThreeDimensionalVector< Eigen::Matrix< double, 1, 3 > >(
        7, numberOfBodies, numberOfBodies, Eigen::Matrix< double, 1, 3 >::Zero( ) );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > > > numericalScalarEihCorrectionsWrtUndergoingPosition = utilities::getThreeDimensionalVector< Eigen::Matrix< double, 1, 3 > >(
        7, numberOfBodies, numberOfBodies, Eigen::Matrix< double, 1, 3 >::Zero( ) );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > > > numericalScalarEihCorrectionsWrtUndergoingVelocity = utilities::getThreeDimensionalVector< Eigen::Matrix< double, 1, 3 > >(
        7, numberOfBodies, numberOfBodies, Eigen::Matrix< double, 1, 3 >::Zero( ) );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 3, 3 > > > > numericalVectorEihCorrectionsWrtExertingPosition = utilities::getThreeDimensionalVector< Eigen::Matrix< double, 3, 3 > >(
        3, numberOfBodies, numberOfBodies, Eigen::Matrix< double, 3, 3 >::Zero( ) );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 3, 3 > > > > numericalVectorEihCorrectionsWrtExertingVelocity = utilities::getThreeDimensionalVector< Eigen::Matrix< double, 3, 3 > >(
        3, numberOfBodies, numberOfBodies, Eigen::Matrix< double, 3, 3 >::Zero( ) );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 3, 3 > > > > numericalVectorEihCorrectionsWrtUndergoingPosition = utilities::getThreeDimensionalVector< Eigen::Matrix< double, 3, 3 > >(
        3, numberOfBodies, numberOfBodies, Eigen::Matrix< double, 3, 3 >::Zero( ) );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 3, 3 > > > > numericalVectorEihCorrectionsWrtUndergoingVelocity = utilities::getThreeDimensionalVector< Eigen::Matrix< double, 3, 3 > >(
        3, numberOfBodies, numberOfBodies, Eigen::Matrix< double, 3, 3 >::Zero( ) );

    std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > >  analyticalTotalPotentialsWrtPosition = eihPartials->getCurrentTotalPotentialWrtPosition( );
    std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > >  analyticalSingleContributionPotentialsWrtPosition = eihPartials->getCurrentLocalPotentialWrtPosition( );
    std::vector< std::vector< Eigen::Matrix< double, 3, 3 > > >  analyticalSinglePointMassAccelerationWrtPosition = eihPartials->getCurrentSinglePointMassAccelerationWrtPosition( );
    std::vector< std::vector< Eigen::Matrix< double, 3, 3 > > >  analyticalTotalPointMassAccelerationWrtPosition = eihPartials->getCurrentTotalPointMassAccelerationWrtPosition( );

    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > > > analyticalScalarEihCorrectionsWrtExertingPosition = utilities::getThreeDimensionalVector< Eigen::Matrix< double, 1, 3 > >(
        7, numberOfBodies, numberOfBodies, Eigen::Matrix< double, 1, 3 >::Zero( ) );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > > > analyticalScalarEihCorrectionsWrtExertingVelocity = utilities::getThreeDimensionalVector< Eigen::Matrix< double, 1, 3 > >(
        7, numberOfBodies, numberOfBodies, Eigen::Matrix< double, 1, 3 >::Zero( ) );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > > > analyticalScalarEihCorrectionsWrtUndergoingPosition = utilities::getThreeDimensionalVector< Eigen::Matrix< double, 1, 3 > >(
        7, numberOfBodies, numberOfBodies, Eigen::Matrix< double, 1, 3 >::Zero( ) );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 1, 3 > > > > analyticalScalarEihCorrectionsWrtUndergoingVelocity = utilities::getThreeDimensionalVector< Eigen::Matrix< double, 1, 3 > >(
        7, numberOfBodies, numberOfBodies, Eigen::Matrix< double, 1, 3 >::Zero( ) );

    std::vector< std::vector< std::vector< Eigen::Matrix< double, 3, 3 > > > > analyticalVectorEihCorrectionsWrtExertingPosition = utilities::getThreeDimensionalVector< Eigen::Matrix< double, 3, 3 > >(
        3, numberOfBodies, numberOfBodies, Eigen::Matrix< double, 3, 3 >::Zero( ) );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 3, 3 > > > > analyticalVectorEihCorrectionsWrtExertingVelocity = utilities::getThreeDimensionalVector< Eigen::Matrix< double, 3, 3 > >(
        3, numberOfBodies, numberOfBodies, Eigen::Matrix< double, 3, 3 >::Zero( ) );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 3, 3 > > > > analyticalVectorEihCorrectionsWrtUndergoingPosition = utilities::getThreeDimensionalVector< Eigen::Matrix< double, 3, 3 > >(
        3, numberOfBodies, numberOfBodies, Eigen::Matrix< double, 3, 3 >::Zero( ) );
    std::vector< std::vector< std::vector< Eigen::Matrix< double, 3, 3 > > > > analyticalVectorEihCorrectionsWrtUndergoingVelocity = utilities::getThreeDimensionalVector< Eigen::Matrix< double, 3, 3 > >(
        3, numberOfBodies, numberOfBodies, Eigen::Matrix< double, 3, 3 >::Zero( ) );

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

                    // Compute analytical d(fs)_{k,ji}/dx_{i}
                    eihPartials->addSingleScalarTermWrtPositionPartial(
                        analyticalScalarEihCorrectionsWrtExertingPosition[ k ][ j ][ i ], j, i, true, k );
                    eihPartials->addSingleScalarTermWrtVelocityPartial(
                        analyticalScalarEihCorrectionsWrtExertingVelocity[ k ][ j ][ i ], j, i, true, k );

                    // Compute analytical d(fs)_{k,ji}/dx_{j}
                    eihPartials->addSingleScalarTermWrtPositionPartial(
                        analyticalScalarEihCorrectionsWrtUndergoingPosition[ k ][ j ][ i ], j, i, false, k );
                    eihPartials->addSingleScalarTermWrtVelocityPartial(
                        analyticalScalarEihCorrectionsWrtUndergoingVelocity[ k ][ j ][ i ], j, i, true, k );

                    if ( k < 3 )
                    {
                        // Compute analytical d(fv)_{k,ji}/dx_{i}
                        eihPartials->addSingleVectorTermWrtPositionPartial(
                            analyticalVectorEihCorrectionsWrtExertingPosition[ k ][ j ][ i ], j, i, true, k );
                        eihPartials->addSingleVectorTermWrtVelocityPartial(
                            analyticalVectorEihCorrectionsWrtExertingVelocity[ k ][ j ][ i ], j, i, true, k );

                        // Compute analytical d(fv)_{k,ji}/dx_{j}
                        eihPartials->addSingleVectorTermWrtPositionPartial(
                            analyticalVectorEihCorrectionsWrtUndergoingPosition[ k ][ j ][ i ], j, i, false, k );
                        eihPartials->addSingleVectorTermWrtVelocityPartial(
                            analyticalVectorEihCorrectionsWrtUndergoingVelocity[ k ][ j ][ i ], j, i, true, k );
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
            numericalTotalPotentialsWrtPosition[ i ][ j ] =
                ( upperturbedTotalPotentials_[ i ][ j ].block( 0, 0, 1, 3 ) - downperturbedTotalPotentials_[ i ][ j ].block( 0, 0, 1, 3 ) ) /
                ( 2.0 * positionPerturbation );


            numericalSingleContributionPotentialsWrtPosition[ i ][ j ] =
                ( upperturbedSingleContrubutionPotentials[ i ][ j ].block( 0, 0, 1, 3 ) - downperturbedSingleContrubutionPotentials[ i ][ j ].block( 0, 0, 1, 3 ) ) /
                ( 2.0 * positionPerturbation );
            numericalSinglePointMassAccelerationWrtPosition[ i ][ j ] =
                ( upperturbedSinglePointmassAcceleration[ i ][ j ].block( 0, 0, 3, 3 ) - downperturbedSinglePointmassAcceleration[ i ][ j ].block( 0, 0, 3, 3 ) ) /
                ( 2.0 * positionPerturbation );
            numericalTotalPointMassAccelerationWrtPosition[ i ][ j ] =
                ( upperturbedTotalPointmassAcceleration[ i ][ j ].block( 0, 0, 3, 3 ) - downperturbedTotalPointmassAcceleration[ i ][ j ].block( 0, 0, 3, 3 ) ) /
                ( 2.0 * positionPerturbation );
            if( i != j )
            {
                for( int k = 0; k < 7; k++ )
                {
                    numericalScalarEihCorrectionsWrtExertingPosition[ k ][ i ][ j ] =
                        ( upperturbedExertingScalarEihCorrections[ k ][ i ][ j ].block( 0, 0, 1, 3 ) - downperturbedExertingScalarEihCorrections[ k ][ i ][ j ].block( 0, 0, 1, 3 ) ) /
                        ( 2.0 * positionPerturbation );

                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        numericalScalarEihCorrectionsWrtExertingPosition[ k ][ i ][ j ].block( 0, 0, 1, 3 ),
                        analyticalScalarEihCorrectionsWrtExertingPosition[ k ][ j ][ i ].block( 0, 0, 1, 3 ),
                        1.0E-4);

                    numericalScalarEihCorrectionsWrtExertingVelocity[ k ][ i ][ j ] =
                        ( upperturbedExertingScalarEihCorrections[ k ][ i ][ j ].block( 0, 3, 1, 3 ) - downperturbedExertingScalarEihCorrections[ k ][ i ][ j ].block( 0, 3, 1, 3 ) ) /
                        ( 2.0 * velocityPerturbation );

                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        numericalScalarEihCorrectionsWrtExertingVelocity[ k ][ i ][ j ].block( 0, 0, 1, 3 ),
                        analyticalScalarEihCorrectionsWrtExertingVelocity[ k ][ j ][ i ].block( 0, 0, 1, 3 ),
                        1.0E-4);

//                    std::cout<<"k "<<k<<" (i,j) "<<i<<" "<<j<<std::endl
//                             <<numericalScalarEihCorrectionsWrtExertingVelocity[ k ][ i ][ j ]<<std::endl
//                             <<analyticalScalarEihCorrectionsWrtExertingVelocity[ k ][ j ][ i ]<<std::endl<<std::endl;

                }
            }

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                numericalTotalPotentialsWrtPosition[ i ][ j ].block( 0, 0, 1, 3 ),
                analyticalTotalPotentialsWrtPosition[ j ][ i ].block( 0, 0, 1, 3 ),
                1.0E-4);

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                numericalTotalPointMassAccelerationWrtPosition[ i ][ j ].block( 0, 0, 3, 3 ),
                analyticalTotalPointMassAccelerationWrtPosition[ j ][ i ].block( 0, 0, 3, 3 ),
                1.0E-4);

            if( i != j )
            {

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    numericalSingleContributionPotentialsWrtPosition[ i ][ j ].block( 0, 0, 1, 3 ),
                    analyticalSingleContributionPotentialsWrtPosition[ j ][ i ].block( 0, 0, 1, 3 ),
                    1.0E-4);

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    numericalSinglePointMassAccelerationWrtPosition[ i ][ j ].block( 0, 0, 3, 3 ),
                    analyticalSinglePointMassAccelerationWrtPosition[ j ][ i ].block( 0, 0, 3, 3 ),
                    1.0E-4);
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
