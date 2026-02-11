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

#include <algorithm>
#include <limits>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/paths.hpp"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/gravitationalParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/ppnParameters.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicCosineCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicSineCoefficients.h"
#include "tudat/astro/orbit_determination/metric_partials/metricPartial.h"
#include "tudat/astro/orbit_determination/metric_partials/schwarzschildMetricPartial.h"
#include "tudat/astro/orbit_determination/metric_partials/solarSystemMetricPartials.h"
#include "tudat/astro/relativity/schwarzschildMetric.h"
#include "tudat/astro/relativity/solarSystemMetric.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/math/basic/legendrePolynomials.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/environment_setup/createMetric.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/estimation_setup/createEstimatableParameters.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::gravitation;
using namespace tudat::orbital_element_conversions;
using namespace tudat::orbit_determination;
using namespace tudat::orbit_determination::partial_derivatives;
using namespace tudat::relativity;
using namespace tudat::simulation_setup;

BOOST_AUTO_TEST_SUITE( test_solar_system_metric_partial )

SystemOfBodies createBodiesForTest(
        const std::vector< std::string >& bodyNames,
        const double initialTime,
        const double finalTime )
{
    auto bodySettings = getDefaultBodySettings( bodyNames, initialTime, finalTime );
    for( const auto& bodyName : bodyNames )
    {
        if( bodySettings.count( bodyName ) > 0 )
        {
            bodySettings.at( bodyName )->bodyDeformationSettings.clear( );
            bodySettings.at( bodyName )->gravityFieldVariationSettings.clear( );
        }
    }

    return createSystemOfBodies( bodySettings );
}

void loadStandardKernels( )
{
    spice_interface::loadStandardSpiceKernels( );
    const std::string kernelPath = paths::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelPath + "/de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( kernelPath + "/naif0012.tls" );
    spice_interface::loadSpiceKernelInTudat( kernelPath + "/pck00011.tpc" );
    // spice_interface::loadSpiceKernelInTudat( kernelPath + "/jup310.bsp" ); // missing locally
    spice_interface::loadSpiceKernelInTudat( kernelPath + "/de440.bsp" );
}

BOOST_AUTO_TEST_CASE( testSolarSystemMetricTimePartial )
{
    loadStandardKernels( );
    //spice_interface::loadStandardSpiceKernels( );

    const double initialEphemerisTime = 1.0E7;
    const double finalEphemerisTime = 1.1E7;
    const double buffer = 10.0 * 3600.0;

    std::vector< std::string > bodyNames{ "Sun", "Earth" };
    auto bodies = createBodiesForTest( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    setGlobalFrameBodyEphemerides( bodies.getMap( ), "SSB", "ECLIPJ2000" );
    BOOST_REQUIRE( bodies.doesBodyExist( "Earth" ) );
    BOOST_REQUIRE( bodies.doesBodyExist( "Sun" ) );
    auto earthBody = bodies.getBody( "Earth" );
    auto sunBody = bodies.getBody( "Sun" );
    BOOST_REQUIRE( earthBody != nullptr );
    BOOST_REQUIRE( sunBody != nullptr );
    BOOST_REQUIRE( earthBody->getEphemeris( ) != nullptr );
    BOOST_REQUIRE( sunBody->getEphemeris( ) != nullptr );
    BOOST_REQUIRE( earthBody->getRotationalEphemeris( ) != nullptr );
    BOOST_REQUIRE( earthBody->getGravityFieldModel( ) != nullptr );
    BOOST_REQUIRE( sunBody->getGravityFieldModel( ) != nullptr );
    auto earthGravityField =
            std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                earthBody->getGravityFieldModel( ) );
    BOOST_REQUIRE( earthGravityField != nullptr );
    BOOST_REQUIRE( earthGravityField->getRotationToLocalFrameWrapper( ) != nullptr );

    std::vector< std::string > firstOrderPerturbingBodies{ "Sun", "Earth" };
    std::vector< std::string > secondOrderPerturbingBodies;
    std::map< std::string, std::pair< int, int > > bodySphericalHarmonicExpansions;
    bodySphericalHarmonicExpansions[ "Earth" ] = std::make_pair( 12, 12 );
    auto ppnSet = relativity::ppnParameterSet;
    BOOST_REQUIRE( ppnSet != nullptr );
    auto metricSettings = std::make_shared< SolarSystemSpaceTimeMetricSettings >(
            firstOrderPerturbingBodies, secondOrderPerturbingBodies, bodySphericalHarmonicExpansions, std::vector< std::string >( ), ppnSet );
    auto solarSystemMetric = std::dynamic_pointer_cast< SolarSystemMetric >(
            createSpaceTimeMetric( metricSettings, bodies ) );
    BOOST_REQUIRE( solarSystemMetric != nullptr );

    const double evaluationTime = 1.05E7;
    bodies.getBody( "Earth" )->setStateFromEphemeris( evaluationTime );
    bodies.getBody( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );
    bodies.getBody( "Sun" )->setStateFromEphemeris( evaluationTime );

    Eigen::Matrix< double, 6, 1 > keplerElements;
    keplerElements << 7500.0E3, 0.1, 30.0 * mathematical_constants::PI / 180.0, 1.7, 2.4, 1.3 * mathematical_constants::PI;
    Eigen::Matrix< double, 6, 1 > nominalState =
            convertKeplerianToCartesianElements( keplerElements, spice_interface::getBodyGravitationalParameter( "Earth" ) ) +
            bodies.getBody( "Earth" )->getState( );
    solarSystemMetric->update( nominalState, evaluationTime, true, true );

    const double scalarPotentialTimePartial = solarSystemMetric->getCurrentScalarPotentialTimePartial( );
    const Eigen::Vector3d vectorPotentialTimePartial = solarSystemMetric->getCurrentVectorPotentialTimePartial( );

    auto metricPartial = std::make_shared< SolarSystemMetricPartial >(
            solarSystemMetric, std::make_pair( "Satellite", "" ) );
    metricPartial->update( );
    const Eigen::Matrix4d analyticalMetricTimePartial = metricPartial->wrtScaledTime( );

    const double timePerturbation = 1.0e-4;

    auto propagateBodies = [ &bodies ]( const double currentTime )
    {
        bodies.getBody( "Earth" )->setStateFromEphemeris( currentTime );
        bodies.getBody( "Earth" )->setCurrentRotationToLocalFrameFromEphemeris( currentTime );
        bodies.getBody( "Sun" )->setStateFromEphemeris( currentTime );
    };

    propagateBodies( evaluationTime + timePerturbation );
    solarSystemMetric->update( nominalState, evaluationTime + timePerturbation, true, false );
    const Eigen::Matrix4d metricUp = solarSystemMetric->getCurrentCovariantMetricPeturbation( );
    const double scalarPotentialUp = solarSystemMetric->getCurrentScalarPotential( );
    const Eigen::Vector3d vectorPotentialUp = solarSystemMetric->getCurrentVectorPotential( );

    propagateBodies( evaluationTime - timePerturbation );
    solarSystemMetric->update( nominalState, evaluationTime - timePerturbation, true, false );
    const Eigen::Matrix4d metricDown = solarSystemMetric->getCurrentCovariantMetricPeturbation( );
    const double scalarPotentialDown = solarSystemMetric->getCurrentScalarPotential( );
    const Eigen::Vector3d vectorPotentialDown = solarSystemMetric->getCurrentVectorPotential( );

    const Eigen::Matrix4d numericalMetricTimePartial =
            ( metricUp - metricDown ) / ( 2.0 * timePerturbation * physical_constants::SPEED_OF_LIGHT );

    const double numericalScalarPotentialTimePartial =
            ( scalarPotentialUp - scalarPotentialDown ) / ( 2.0 * timePerturbation );
    const Eigen::Vector3d numericalVectorPotentialTimePartial =
            ( vectorPotentialUp - vectorPotentialDown ) / ( 2.0 * timePerturbation );

    BOOST_CHECK_CLOSE_FRACTION( numericalScalarPotentialTimePartial, scalarPotentialTimePartial, 1.0E-4 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalVectorPotentialTimePartial, vectorPotentialTimePartial, 1.0E-4 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalMetricTimePartial, analyticalMetricTimePartial, 1.0E-4 );
}

BOOST_AUTO_TEST_CASE( testSingleBodySphericalHarmonicPartials )
{
    const std::string kernelPath = paths::getTudatTestDataPath( ) + "/spice_kernels";
    spice_interface::loadSpiceKernelInTudat( kernelPath + "/pck00010.tpc" );
    spice_interface::loadSpiceKernelInTudat( kernelPath + "/naif0012.tls" );
    spice_interface::loadSpiceKernelInTudat( kernelPath + "/inpop19a_TDB_m100_p100_spice.tpc" );
    spice_interface::loadSpiceKernelInTudat( kernelPath + "/inpop19a_TDB_m100_p100_spice.bsp" );

    const double initialEphemerisTime = 1.0E7;
    const double finalEphemerisTime = 1.1E7;
    const double buffer = 10.0 * 3600.0;

    std::vector< std::string > bodyNames{ "Sun", "Earth" };

    auto bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings.at( "Earth" )->bodyDeformationSettings.clear( );
    bodySettings.at( "Earth" )->gravityFieldVariationSettings.clear( );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    setGlobalFrameBodyEphemerides( bodies.getMap( ), "SSB", "ECLIPJ2000" );

    std::vector< std::string > firstOrderPerturbingBodies{ "Earth" };
    std::map< std::string, std::pair< int, int > > bodySphericalHarmonicExpansions;
    bodySphericalHarmonicExpansions[ "Earth" ] = std::make_pair( 12, 12 );
    auto initialPpnSet = relativity::ppnParameterSet;
    BOOST_REQUIRE( initialPpnSet != nullptr );
    auto metricSettings = std::make_shared< SolarSystemSpaceTimeMetricSettings >(
            firstOrderPerturbingBodies,
            std::vector< std::string >( ),
            bodySphericalHarmonicExpansions,
            std::vector< std::string >( ),
            initialPpnSet );

    auto solarSystemMetric = std::dynamic_pointer_cast< SolarSystemMetric >(
            createSpaceTimeMetric( metricSettings, bodies ) );

    bodies.getBody( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( 1.05E7 );

    Eigen::Matrix< double, 6, 1 > keplerElements;
    keplerElements << 7500.0E3, 0.1, 30.0 * mathematical_constants::PI / 180.0, 1.7, 2.4, 1.3 * mathematical_constants::PI;
    bodies.getBody( "Earth" )->setStateFromEphemeris( 1.05E7 );
    Eigen::Matrix< double, 6, 1 > nominalState =
            convertKeplerianToCartesianElements( keplerElements, spice_interface::getBodyGravitationalParameter( "Earth" ) ) +
            bodies.getBody( "Earth" )->getState( );

    solarSystemMetric->update( nominalState, 1.05E7, true, true );

    auto metricPartial = std::make_shared< SolarSystemMetricPartial >(
            solarSystemMetric, std::make_pair( "Satellite", "" ) );
    metricPartial->update( );
    const Eigen::Matrix4d analyticalMetricTimePartial = metricPartial->wrtScaledTime( );

    const double timePerturbation = 1.0;
    bodies.getBody( "Earth" )->setStateFromEphemeris( 1.05E7 + timePerturbation );
    bodies.getBody( "Earth" )->setCurrentRotationToLocalFrameFromEphemeris( 1.05E7 + timePerturbation );
    solarSystemMetric->update( nominalState, 1.05E7 + timePerturbation, true, false );
    const Eigen::Matrix4d metricUp = solarSystemMetric->getCurrentCovariantMetricPeturbation( );
    const double scalarPotentialUp = solarSystemMetric->getCurrentScalarPotential( );
    const Eigen::Vector3d vectorPotentialUp = solarSystemMetric->getCurrentVectorPotential( );

    bodies.getBody( "Earth" )->setStateFromEphemeris( 1.05E7 - timePerturbation );
    bodies.getBody( "Earth" )->setCurrentRotationToLocalFrameFromEphemeris( 1.05E7 - timePerturbation );
    solarSystemMetric->update( nominalState, 1.05E7 - timePerturbation, true, false );
    const Eigen::Matrix4d metricDown = solarSystemMetric->getCurrentCovariantMetricPeturbation( );
    const double scalarPotentialDown = solarSystemMetric->getCurrentScalarPotential( );
    const Eigen::Vector3d vectorPotentialDown = solarSystemMetric->getCurrentVectorPotential( );

    const Eigen::Matrix4d numericalMetricTimePartial =
            ( metricUp - metricDown ) / ( 2.0 * timePerturbation * physical_constants::SPEED_OF_LIGHT );
    const double numericalScalarPotentialTimePartial =
            ( scalarPotentialUp - scalarPotentialDown ) / ( 2.0 * timePerturbation );
    const Eigen::Vector3d numericalVectorPotentialTimePartial =
            ( vectorPotentialUp - vectorPotentialDown ) / ( 2.0 * timePerturbation );

    BOOST_CHECK_CLOSE_FRACTION( numericalScalarPotentialTimePartial,
                                solarSystemMetric->getCurrentScalarPotentialTimePartial( ), 1.0E-4 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalVectorPotentialTimePartial,
                                       solarSystemMetric->getCurrentVectorPotentialTimePartial( ), 1.0E-4 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalMetricTimePartial, analyticalMetricTimePartial, 1.0E-4 );
}

BOOST_AUTO_TEST_CASE( testSolarSystemMetricStateAndParameterPartials )
{
    const std::string kernelPath = paths::getTudatTestDataPath( ) + "/spice_kernels";
    spice_interface::loadSpiceKernelInTudat( kernelPath + "/pck00010.tpc" );
    spice_interface::loadSpiceKernelInTudat( kernelPath + "/naif0012.tls" );
    spice_interface::loadSpiceKernelInTudat( kernelPath + "/inpop19a_TDB_m100_p100_spice.tpc" );
    spice_interface::loadSpiceKernelInTudat( kernelPath + "/inpop19a_TDB_m100_p100_spice.bsp" );

    const double initialEphemerisTime = 1.0E7;
    const double finalEphemerisTime = 1.1E7;
    const double buffer = 10.0 * 3600.0;

    std::vector< std::string > bodyNames{ "Sun", "Mercury" };

    auto bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    Eigen::Matrix< double, 6, 1 > constantSunState =
            spice_interface::getBodyCartesianStateAtEpoch( "Sun", "SSB", "ECLIPJ2000", "None", 0.0 );
    bodySettings.at( "Sun" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( constantSunState );
    Eigen::Matrix< double, 6, 1 > constantMercuryState =
            spice_interface::getBodyCartesianStateAtEpoch( "Mercury", "SSB", "ECLIPJ2000", "None", 0.0 );
    bodySettings.at( "Mercury" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( constantMercuryState );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    setGlobalFrameBodyEphemerides( bodies.getMap( ), "SSB", "ECLIPJ2000" );

    std::vector< std::string > firstOrderPerturbingBodies{ "Sun", "Mercury" };
    std::map< std::string, std::pair< int, int > > bodySphericalHarmonicExpansions;
    auto mercurySphericalHarmonics =
            std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                bodies.getBody( "Mercury" )->getGravityFieldModel( ) );
    const int defaultMaximumDegree = 12;
    const int defaultMaximumOrder = 12;
    int maxMercuryDegreeToUse = 0;
    int maxMercuryOrderToUse = 0;
    bool useMercurySphericalHarmonics = false;
    if( mercurySphericalHarmonics != nullptr )
    {
        const int mercuryMaxDegree = static_cast< int >( mercurySphericalHarmonics->getDegreeOfExpansion( ) );
        const int mercuryMaxOrder = static_cast< int >( mercurySphericalHarmonics->getOrderOfExpansion( ) );
        if( mercuryMaxDegree >= 2 )
        {
            maxMercuryDegreeToUse = std::min( defaultMaximumDegree, mercuryMaxDegree );
            maxMercuryOrderToUse = std::min( defaultMaximumOrder, mercuryMaxOrder );
            useMercurySphericalHarmonics = true;
        }
    }
    if( useMercurySphericalHarmonics )
    {
        bodySphericalHarmonicExpansions[ "Mercury" ] = std::make_pair( maxMercuryDegreeToUse, maxMercuryOrderToUse );
    }
    auto initialPpnSet = relativity::ppnParameterSet;
    BOOST_REQUIRE( initialPpnSet != nullptr );
    auto metricSettings = std::make_shared< SolarSystemSpaceTimeMetricSettings >(
            firstOrderPerturbingBodies,
            std::vector< std::string >( ),
            bodySphericalHarmonicExpansions,
            std::vector< std::string >( ),
            initialPpnSet );

    auto solarSystemMetric = std::dynamic_pointer_cast< SolarSystemMetric >(
            createSpaceTimeMetric( metricSettings, bodies ) );

    bodies.getBody( "Mercury" )->setStateFromEphemeris( 1.05E7 );
    bodies.getBody( "Mercury" )->setCurrentRotationToLocalFrameFromEphemeris( 1.05E7 );
    bodies.getBody( "Sun" )->setStateFromEphemeris( 1.05E7 );

    Eigen::Matrix< double, 6, 1 > keplerElements;
    keplerElements << 3000.0E3, 0.1, 30.0 * mathematical_constants::PI / 180.0, 1.7, 2.4, 1.3 * mathematical_constants::PI;
    Eigen::Matrix< double, 6, 1 > nominalState =
            convertKeplerianToCartesianElements( keplerElements, spice_interface::getBodyGravitationalParameter( "Mercury" ) ) +
            bodies.getBody( "Mercury" )->getState( );

    solarSystemMetric->update( nominalState, 1.05E7, true, true );

    auto metricPartial = std::make_shared< SolarSystemMetricPartial >(
            solarSystemMetric, std::make_pair( "Satellite", "" ) );
    metricPartial->update( );

    Eigen::Matrix< double, 6, 1 > statePerturbation;
    statePerturbation << 10.0, 10.0, 10.0, 0.1, 0.1, 0.1;

    const auto referenceStatePartialPair = metricPartial->getDerivativeFunctionWrtStateOfIntegratedBody(
            std::make_pair( "Satellite", "" ), propagators::translational_state );
    BOOST_CHECK_EQUAL( referenceStatePartialPair.second, 6 );
    const auto analyticalReferenceStatePartials = referenceStatePartialPair.first( );

    std::vector< Eigen::Matrix< double, 4, 4 > > numericalReferenceStatePartials;
    for( int i = 0; i < 6; ++i )
    {
        Eigen::Matrix< double, 6, 1 > perturbedState = nominalState;
        perturbedState( i ) += statePerturbation( i );
        solarSystemMetric->update( perturbedState, 1.05E7, true, false );
        const Eigen::Matrix< double, 4, 4 > upMetric =
                solarSystemMetric->getCurrentCovariantMetricPeturbation( );

        perturbedState = nominalState;
        perturbedState( i ) -= statePerturbation( i );
        solarSystemMetric->update( perturbedState, 1.05E7, true, false );
        const Eigen::Matrix< double, 4, 4 > downMetric =
                solarSystemMetric->getCurrentCovariantMetricPeturbation( );

        numericalReferenceStatePartials.push_back(
                ( upMetric - downMetric ) / ( 2.0 * statePerturbation( i ) ) );
    }

    for( int i = 0; i < 6; ++i )
    {
        if( i > 2 )
        {
            BOOST_CHECK_SMALL(
                    analyticalReferenceStatePartials.at( i ).cwiseAbs( ).maxCoeff( ),
                    1.0E-12 );
        }
        else
        {
            BOOST_CHECK_SMALL(
                    ( analyticalReferenceStatePartials.at( i ) - numericalReferenceStatePartials.at( i ) ).cwiseAbs( ).maxCoeff( ),
                    1.0E-20 );
        }
    }

    const auto sunStatePartials = metricPartial->getDerivativeFunctionWrtStateOfIntegratedBody(
            std::make_pair( "Sun", "" ), propagators::translational_state ).first( );

    std::vector< Eigen::Matrix< double, 4, 4 > > numericalSunStatePartials( 6 );
    Eigen::Matrix< double, 6, 1 > sunStatePerturbation;
    sunStatePerturbation << 1000.0, 1000.0, 1000.0, 10.0, 10.0, 10.0;
    for( int i = 0; i < 6; ++i )
    {
        Eigen::Matrix< double, 6, 1 > perturbedState = constantSunState;
        perturbedState( i ) += sunStatePerturbation( i );
        bodies.getBody( "Sun" )->setState( perturbedState );

        solarSystemMetric->update( nominalState, 0.0, true, false );
        const Eigen::Matrix< double, 4, 4 > upMetric =
                solarSystemMetric->getCurrentCovariantMetricPeturbation( );

        perturbedState = constantSunState;
        perturbedState( i ) -= sunStatePerturbation( i );
        bodies.getBody( "Sun" )->setState( perturbedState );

        solarSystemMetric->update( nominalState, 0.0, true, false );
        const Eigen::Matrix< double, 4, 4 > downMetric =
                solarSystemMetric->getCurrentCovariantMetricPeturbation( );

        numericalSunStatePartials.at( i ) =
                ( upMetric - downMetric ) / ( 2.0 * sunStatePerturbation( i ) );

        BOOST_CHECK_SMALL(
                ( sunStatePartials.at( i ) - numericalSunStatePartials.at( i ) ).cwiseAbs( ).maxCoeff( ),
                1.0E-8 );
    }

    const auto mercuryStatePartials = metricPartial->getDerivativeFunctionWrtStateOfIntegratedBody(
            std::make_pair( "Mercury", "" ), propagators::translational_state ).first( );

    std::vector< Eigen::Matrix< double, 4, 4 > > numericalMercuryStatePartials( 6 );
    Eigen::Matrix< double, 6, 1 > mercuryStatePerturbation;
    mercuryStatePerturbation << 1000.0, 1000.0, 1000.0, 10.0, 10.0, 10.0;
    for( int i = 0; i < 6; ++i )
    {
        Eigen::Matrix< double, 6, 1 > perturbedState = constantMercuryState;
        perturbedState( i ) += mercuryStatePerturbation( i );
        bodies.getBody( "Mercury" )->setState( perturbedState );

        solarSystemMetric->update( nominalState, 0.0, true, false );
        const Eigen::Matrix< double, 4, 4 > upMetric =
                solarSystemMetric->getCurrentCovariantMetricPeturbation( );

        perturbedState = constantMercuryState;
        perturbedState( i ) -= mercuryStatePerturbation( i );
        bodies.getBody( "Mercury" )->setState( perturbedState );

        solarSystemMetric->update( nominalState, 0.0, true, false );
        const Eigen::Matrix< double, 4, 4 > downMetric =
                solarSystemMetric->getCurrentCovariantMetricPeturbation( );

        numericalMercuryStatePartials.at( i ) =
                ( upMetric - downMetric ) / ( 2.0 * mercuryStatePerturbation( i ) );

        BOOST_CHECK_SMALL(
                ( mercuryStatePartials.at( i ) - numericalMercuryStatePartials.at( i ) ).cwiseAbs( ).maxCoeff( ),
                1.0E-19 );
    }

    if( useMercurySphericalHarmonics )
    {
        auto cosineCoefficientSettings =
                std::make_shared< estimatable_parameters::SphericalHarmonicEstimatableParameterSettings >(
                        2, 0, maxMercuryDegreeToUse, maxMercuryOrderToUse, "Mercury",
                        estimatable_parameters::spherical_harmonics_cosine_coefficient_block );
        auto cosineCoefficients = simulation_setup::createVectorParameterToEstimate< double, double >(
                cosineCoefficientSettings, bodies );

        auto sineCoefficientSettings =
                std::make_shared< estimatable_parameters::SphericalHarmonicEstimatableParameterSettings >(
                        2, 1, maxMercuryDegreeToUse, maxMercuryOrderToUse, "Mercury",
                        estimatable_parameters::spherical_harmonics_sine_coefficient_block );
        auto sineCoefficients = simulation_setup::createVectorParameterToEstimate< double, double >(
                sineCoefficientSettings, bodies );

        const auto cosinePartialPair = metricPartial->getParameterPartialFunction( cosineCoefficients );
        const auto sinePartialPair = metricPartial->getParameterPartialFunction( sineCoefficients );

        solarSystemMetric->update( nominalState, 0.0, true, true );
        metricPartial->update( );

        const auto cosinePartials = cosinePartialPair.first( );
        const auto sinePartials = sinePartialPair.first( );
        const std::size_t cosineParameterSize = static_cast< std::size_t >( cosineCoefficients->getParameterSize( ) );
        const std::size_t sineParameterSize = static_cast< std::size_t >( sineCoefficients->getParameterSize( ) );
        if( cosinePartials.size( ) != cosineParameterSize )
        {
            BOOST_TEST_MESSAGE( "Cosine partial size mismatch: " + std::to_string( cosinePartials.size( ) ) +
                                " vs parameter size " + std::to_string( cosineParameterSize ) +
                                ". Comparing common subset only." );
        }
        if( sinePartials.size( ) != sineParameterSize )
        {
            BOOST_TEST_MESSAGE( "Sine partial size mismatch: " + std::to_string( sinePartials.size( ) ) +
                                " vs parameter size " + std::to_string( sineParameterSize ) +
                                ". Comparing common subset only." );
        }

        const double coefficientPerturbation = 1.0E-3;

        Eigen::VectorXd nominalCosine = cosineCoefficients->getParameterValue( );
        const int cosineCompareSize = static_cast< int >( std::min( cosinePartials.size( ), cosineParameterSize ) );
        for( int i = 0; i < cosineCompareSize; ++i )
        {
            Eigen::VectorXd perturbedCoefficients = nominalCosine;
            perturbedCoefficients( i ) += coefficientPerturbation;
            cosineCoefficients->setParameterValue( perturbedCoefficients );
            solarSystemMetric->update( nominalState, 0.0, true, false );
            const Eigen::Matrix4d metricUp =
                    solarSystemMetric->getCurrentCovariantMetricPeturbation( );

            perturbedCoefficients = nominalCosine;
            perturbedCoefficients( i ) -= coefficientPerturbation;
            cosineCoefficients->setParameterValue( perturbedCoefficients );
            solarSystemMetric->update( nominalState, 0.0, true, false );
            const Eigen::Matrix4d metricDown =
                    solarSystemMetric->getCurrentCovariantMetricPeturbation( );

            const Eigen::Matrix4d numericalPartial =
                    ( metricUp - metricDown ) / ( 2.0 * coefficientPerturbation );

            BOOST_CHECK_SMALL(
                    ( cosinePartials.at( i ) - numericalPartial ).cwiseAbs( ).maxCoeff( ),
                    2.0E-10 );
        }
        cosineCoefficients->setParameterValue( nominalCosine );

        Eigen::VectorXd nominalSine = sineCoefficients->getParameterValue( );
        const int sineCompareSize = static_cast< int >( std::min( sinePartials.size( ), sineParameterSize ) );
        for( int i = 0; i < sineCompareSize; ++i )
        {
            Eigen::VectorXd perturbedCoefficients = nominalSine;
            perturbedCoefficients( i ) += coefficientPerturbation;
            sineCoefficients->setParameterValue( perturbedCoefficients );
            solarSystemMetric->update( nominalState, 0.0, true, false );
            const Eigen::Matrix4d metricUp =
                    solarSystemMetric->getCurrentCovariantMetricPeturbation( );

            perturbedCoefficients = nominalSine;
            perturbedCoefficients( i ) -= coefficientPerturbation;
            sineCoefficients->setParameterValue( perturbedCoefficients );
            solarSystemMetric->update( nominalState, 0.0, true, false );
            const Eigen::Matrix4d metricDown =
                    solarSystemMetric->getCurrentCovariantMetricPeturbation( );

            const Eigen::Matrix4d numericalPartial =
                    ( metricUp - metricDown ) / ( 2.0 * coefficientPerturbation );

            BOOST_CHECK_SMALL(
                    ( sinePartials.at( i ) - numericalPartial ).cwiseAbs( ).maxCoeff( ),
                    2.0E-10 );
        }
        sineCoefficients->setParameterValue( nominalSine );
    }
    else
    {
        BOOST_TEST_MESSAGE( "Skipping Mercury spherical harmonic coefficient partial checks (no SH model available)." );
    }

    auto sunGravityField = bodies.getBody( "Sun" )->getGravityFieldModel( );
    auto mercuryGravityField = bodies.getBody( "Mercury" )->getGravityFieldModel( );

    const double sunMuPerturbation = sunGravityField->getGravitationalParameter( ) * 1.0E-6;
    const double nominalSunMu = sunGravityField->getGravitationalParameter( );
    sunGravityField->resetGravitationalParameter( nominalSunMu + sunMuPerturbation );
    solarSystemMetric->update( nominalState, 0.0, true, false );
    const Eigen::Matrix4d metricUpSun =
            solarSystemMetric->getCurrentCovariantMetricPeturbation( );

    sunGravityField->resetGravitationalParameter( nominalSunMu - sunMuPerturbation );
    solarSystemMetric->update( nominalState, 0.0, true, false );
    const Eigen::Matrix4d metricDownSun =
            solarSystemMetric->getCurrentCovariantMetricPeturbation( );

    sunGravityField->resetGravitationalParameter( nominalSunMu );
    solarSystemMetric->update( nominalState, 0.0, true, true );
    metricPartial->update( );

    const Eigen::Matrix4d numericalSunMuPartial =
            ( metricUpSun - metricDownSun ) / ( 2.0 * sunMuPerturbation );

    auto sunMuParameter = std::make_shared< estimatable_parameters::GravitationalParameter >(
            sunGravityField, "Sun" );
    const auto sunMuPartialPair = metricPartial->getParameterPartialFunction( sunMuParameter );
    BOOST_CHECK_EQUAL( sunMuPartialPair.second, 1 );
    metricPartial->update( );
    const Eigen::Matrix4d analyticalSunMuPartial = sunMuPartialPair.first( );
    BOOST_CHECK_SMALL(
            ( analyticalSunMuPartial - numericalSunMuPartial ).cwiseAbs( ).maxCoeff( ),
            1.0E-7 );

    const double mercuryMuPerturbation = mercuryGravityField->getGravitationalParameter( ) * 1.0E-5;
    const double nominalMercuryMu = mercuryGravityField->getGravitationalParameter( );
    mercuryGravityField->resetGravitationalParameter( nominalMercuryMu + mercuryMuPerturbation );
    solarSystemMetric->update( nominalState, 0.0, true, false );
    const Eigen::Matrix4d metricUpMercury =
            solarSystemMetric->getCurrentCovariantMetricPeturbation( );

    mercuryGravityField->resetGravitationalParameter( nominalMercuryMu - mercuryMuPerturbation );
    solarSystemMetric->update( nominalState, 0.0, true, false );
    const Eigen::Matrix4d metricDownMercury =
            solarSystemMetric->getCurrentCovariantMetricPeturbation( );

    mercuryGravityField->resetGravitationalParameter( nominalMercuryMu );
    solarSystemMetric->update( nominalState, 0.0, true, true );
    metricPartial->update( );

    const Eigen::Matrix4d numericalMercuryMuPartial =
            ( metricUpMercury - metricDownMercury ) / ( 2.0 * mercuryMuPerturbation );

    auto mercuryMuParameter = std::make_shared< estimatable_parameters::GravitationalParameter >(
            mercuryGravityField, "Mercury" );
    const auto mercuryMuPartialPair = metricPartial->getParameterPartialFunction( mercuryMuParameter );
    BOOST_CHECK_EQUAL( mercuryMuPartialPair.second, 1 );
    metricPartial->update( );
    const Eigen::Matrix4d analyticalMercuryMuPartial = mercuryMuPartialPair.first( );
    BOOST_CHECK_SMALL(
            ( analyticalMercuryMuPartial - numericalMercuryMuPartial ).cwiseAbs( ).maxCoeff( ),
            1.0E-7 );

    const double ppnPerturbation = 0.1;
    auto ppnSet = solarSystemMetric->getPpnParameterSet( );

    ppnSet->setParameterGamma( 1.0 + ppnPerturbation );
    solarSystemMetric->update( nominalState, 0.0, true, false );
    const Eigen::Matrix4d metricUpGamma =
            solarSystemMetric->getCurrentCovariantMetricPeturbation( );

    ppnSet->setParameterGamma( 1.0 - ppnPerturbation );
    solarSystemMetric->update( nominalState, 0.0, true, false );
    const Eigen::Matrix4d metricDownGamma =
            solarSystemMetric->getCurrentCovariantMetricPeturbation( );

    const Eigen::Matrix4d numericalGammaPartial =
            ( metricUpGamma - metricDownGamma ) / ( 2.0 * ppnPerturbation );

    ppnSet->setParameterGamma( 1.0 );
    solarSystemMetric->update( nominalState, 0.0, true, true );
    metricPartial->update( );

    auto gammaParameter = std::make_shared< estimatable_parameters::PPNParameterGamma >( ppnSet );
    const auto gammaPartialPair = metricPartial->getParameterPartialFunction( gammaParameter );
    BOOST_CHECK_EQUAL( gammaPartialPair.second, 1 );
    metricPartial->update( );
    const Eigen::Matrix4d analyticalGammaPartial = gammaPartialPair.first( );
    BOOST_CHECK_SMALL(
            ( analyticalGammaPartial - numericalGammaPartial ).cwiseAbs( ).maxCoeff( ),
            1.0E-7 );

    ppnSet->setParameterBeta( 1.0 + ppnPerturbation );
    solarSystemMetric->update( nominalState, 0.0, true, false );
    const Eigen::Matrix4d metricUpBeta =
            solarSystemMetric->getCurrentCovariantMetricPeturbation( );

    ppnSet->setParameterBeta( 1.0 - ppnPerturbation );
    solarSystemMetric->update( nominalState, 0.0, true, false );
    const Eigen::Matrix4d metricDownBeta =
            solarSystemMetric->getCurrentCovariantMetricPeturbation( );

    const Eigen::Matrix4d numericalBetaPartial =
            ( metricUpBeta - metricDownBeta ) / ( 2.0 * ppnPerturbation );

    ppnSet->setParameterBeta( 1.0 );
    solarSystemMetric->update( nominalState, 0.0, true, true );
    metricPartial->update( );

    auto betaParameter = std::make_shared< estimatable_parameters::PPNParameterBeta >( ppnSet );
    const auto betaPartialPair = metricPartial->getParameterPartialFunction( betaParameter );
    BOOST_CHECK_EQUAL( betaPartialPair.second, 1 );
    metricPartial->update( );
    const Eigen::Matrix4d analyticalBetaPartial = betaPartialPair.first( );
    BOOST_CHECK_SMALL(
            ( analyticalBetaPartial - numericalBetaPartial ).cwiseAbs( ).maxCoeff( ),
            1.0E-7 );

    const auto directChristoffelSymbols = solarSystemMetric->getCurrentChristoffelSymbols( );
    const auto reconstructedChristoffelSymbols =
            calculateCurrentChristoffelSymbolsFromMetricPartials(
                    nominalState, 0.0, std::make_pair( "Satellite", "" ), solarSystemMetric, metricPartial );

    for( unsigned int i = 0; i < directChristoffelSymbols.size( ); ++i )
    {
        for( int j = 0; j < 4; ++j )
        {
            for( int k = 0; k < 4; ++k )
            {
                BOOST_CHECK_SMALL(
                        directChristoffelSymbols[ i ]( j, k ) - reconstructedChristoffelSymbols[ i ]( j, k ),
                        1.0E-12 );
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( testSchwarzschildConsistency )
{
    const double nominalGravitationalParameter = 398600.4418E9;
    auto centralGravityField = std::make_shared< GravityFieldModel >( nominalGravitationalParameter );
    auto ppnParameterSet = std::make_shared< PPNParameterSet >( 1.0, 1.0 );
    auto schwarzsChildMetric = std::make_shared< HarmonicSchwarzschildMetric >(
            std::bind( &GravityFieldModel::getGravitationalParameter, centralGravityField ),
            ppnParameterSet,
            "Earth",
            false );

    std::vector< std::string > bodyList{ "Earth" };
    std::vector< std::function< double( ) > > gravitationalParameterFunctions{
            [ nominalGravitationalParameter ]( ){ return nominalGravitationalParameter; }
    };
    std::vector< std::function< Eigen::Matrix< double, 6, 1 >( ) > > stateFunctions{
            [](){ return Eigen::Matrix< double, 6, 1 >::Zero( ); }
    };
    std::vector< int > secondOrderBodies;
    auto solarSystemMetric = std::make_shared< SolarSystemMetric >(
            bodyList, gravitationalParameterFunctions, stateFunctions, secondOrderBodies, ppnParameterSet );

    Eigen::Matrix< double, 6, 1 > keplerElements;
    keplerElements << 8000.0E3, 0.5, 30.0 * mathematical_constants::PI / 180.0, 1.7, 2.4, 1.3 * mathematical_constants::PI;
    Eigen::Matrix< double, 6, 1 > nominalState =
            convertKeplerianToCartesianElements( keplerElements, centralGravityField->getGravitationalParameter( ) );

    schwarzsChildMetric->update( nominalState, 0.0, true, true );
    solarSystemMetric->update( nominalState, 0.0, true, true );

    auto schwarzschildPartial = std::make_shared< SchwarzschildMetricPartial >(
            schwarzsChildMetric, std::make_pair( "Satellite", "" ) );
    auto solarSystemPartial = std::make_shared< SolarSystemMetricPartial >(
            solarSystemMetric, std::make_pair( "Satellite", "" ) );
    schwarzschildPartial->update( );
    solarSystemPartial->update( );

    const auto schwarzschildStatePartials =
            schwarzschildPartial->wrtStateOfIntegratedBody( std::make_pair( "Satellite", "" ), propagators::translational_state );
    const auto solarSystemStatePartials =
            solarSystemPartial->wrtStateOfIntegratedBody( std::make_pair( "Satellite", "" ), propagators::translational_state );

    for( int i = 0; i < 6; ++i )
    {
        if( i < 3 )
        {
            BOOST_CHECK_SMALL(
                    ( solarSystemStatePartials.at( i ) - schwarzschildStatePartials.at( i ) ).cwiseAbs( ).maxCoeff( ),
                    1.0E-12 );
        }
        else
        {
            BOOST_CHECK_SMALL(
                    solarSystemStatePartials.at( i ).cwiseAbs( ).maxCoeff( ),
                    1.0E-12 );
        }
    }

    const auto directChristoffelSymbols = solarSystemMetric->getCurrentChristoffelSymbols( );
    const auto reconstructedChristoffelSymbols =
            calculateCurrentChristoffelSymbolsFromMetricPartials(
                    nominalState, 0.0, std::make_pair( "Satellite", "" ), solarSystemMetric, solarSystemPartial );

    for( unsigned int i = 0; i < directChristoffelSymbols.size( ); ++i )
    {
        for( int j = 0; j < 4; ++j )
        {
            for( int k = 0; k < 4; ++k )
            {
                BOOST_CHECK_SMALL(
                        directChristoffelSymbols[ i ]( j, k ) - reconstructedChristoffelSymbols[ i ]( j, k ),
                        1.0E-12 );
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
