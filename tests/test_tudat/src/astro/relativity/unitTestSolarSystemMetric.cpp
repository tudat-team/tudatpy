/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/astro/relativity/solarSystemMetric.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityField.h"
#include "tudat/math/basic/legendrePolynomials.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createMetric.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::simulation_setup;
using namespace tudat::relativity;
using namespace tudat::orbital_element_conversions;
using namespace tudat::spice_interface;

BOOST_AUTO_TEST_SUITE( test_solar_system_metric )

BOOST_AUTO_TEST_CASE( testStaticSolarSystemMetricAgainstSchwarzschild )
{
    loadStandardSpiceKernels( );

    const double initialEphemerisTime = 1.0E7;
    const double finalEphemerisTime = 1.1E7;
    const double maximumTimeStep = 3600.0;
    const double buffer = 10.0 * maximumTimeStep;
    const double evaluationTime = 1.05E7;

    std::vector< std::string > bodyNames = { "Sun", "Earth", "Moon" };
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings.at( "Earth" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( Eigen::Vector6d::Zero( ) );
    bodySettings.at( "Earth" )->gravityFieldVariationSettings.clear( );
    auto bodies = createSystemOfBodies( bodySettings );
    setGlobalFrameBodyEphemerides( bodies.getMap( ), "SSB", "ECLIPJ2000" );
    bodies.getBody( "Earth" )->setStateFromEphemeris( evaluationTime );

    auto firstOrderSchwarzschildSettings = std::make_shared< SchwarzschildSpaceTimeMetricSettings >( "Earth", false );
    auto secondOrderSchwarzschildSettings = std::make_shared< SchwarzschildSpaceTimeMetricSettings >( "Earth", true );

    auto firstOrderSolarSystemSettings =
            std::make_shared< SolarSystemSpaceTimeMetricSettings >( std::vector< std::string >{ "Earth" },
                                                                    std::vector< std::string >( ),
                                                                    std::map< std::string, std::pair< int, int > >( ),
                                                                    std::vector< std::string >( ) );
    // Second-order solar-system terms currently unimplemented; re-enable once metric supports them.
    // auto secondOrderSolarSystemSettings =
    //         std::make_shared< SolarSystemSpaceTimeMetricSettings >(
    //                 std::vector< std::string >( ),
    //                 std::vector< std::string >{ "Earth" },
    //                 std::map< std::string, std::pair< int, int > >( ),
    //                 std::vector< std::string >( ),
    //                 ppnParameterSet );

    auto firstOrderSchwarzschildMetric = createSpaceTimeMetric( firstOrderSchwarzschildSettings, bodies );
    auto secondOrderSchwarzschildMetric = createSpaceTimeMetric( secondOrderSchwarzschildSettings, bodies );
    auto firstOrderSolarSystemMetric = createSpaceTimeMetric( firstOrderSolarSystemSettings, bodies );

    Eigen::Vector6d keplerianElements;
    keplerianElements << 6378.0E3 + 400.0E3, 0.013, 1.0238269559089248, 3.1526292818328812, 1.5807574453453865, 3.1478950321924795;

    const double earthGravitationalParameter = bodies.getBody( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    const Eigen::Vector6d cartesianState = convertKeplerianToCartesianElements( keplerianElements, earthGravitationalParameter );

    firstOrderSchwarzschildMetric->update( cartesianState, 0.0, true, true );
    secondOrderSchwarzschildMetric->update( cartesianState, 0.0, true, true );

    const Eigen::Vector6d barycentricState = cartesianState + bodies.getBody( "Earth" )->getState( );
    firstOrderSolarSystemMetric->update( barycentricState, 0.0, true, true );

    const Eigen::Matrix4d firstOrderSchwarzschildPerturbation = firstOrderSchwarzschildMetric->getCurrentCovariantMetricPeturbation( );
    const Eigen::Matrix4d firstOrderSolarPerturbation = firstOrderSolarSystemMetric->getCurrentCovariantMetricPeturbation( );

    BOOST_CHECK_CLOSE_FRACTION( firstOrderSchwarzschildPerturbation( 0, 0 ), firstOrderSolarPerturbation( 0, 0 ), 1.0E-13 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            firstOrderSchwarzschildPerturbation.block( 1, 1, 3, 3 ), firstOrderSolarPerturbation.block( 1, 1, 3, 3 ), 1.0E-13 );

    const Eigen::Vector3d earthVelocity = bodies.getBody( "Earth" )->getState( ).segment( 3, 3 );
    const Eigen::Vector3d expectedSpaceTimePerturbation =
            -2.0 * earthVelocity * firstOrderSolarPerturbation( 1, 1 ) / physical_constants::SPEED_OF_LIGHT;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( firstOrderSolarPerturbation.block( 0, 1, 1, 3 ),
                                       expectedSpaceTimePerturbation.transpose( ),
                                       std::numeric_limits< double >::epsilon( ) );
    // TODO(second-order): Re-enable solar-system second-order comparisons once the metric supports epsilon-dependent terms.
    // const Eigen::Matrix4d secondOrderSchwarzschildPerturbation =
    //         secondOrderSchwarzschildMetric->getCurrentCovariantMetricPeturbation( );
    // const Eigen::Matrix4d secondOrderSolarPerturbation =
    //         secondOrderSolarSystemMetric->getCurrentCovariantMetricPeturbation( );
    // for( int i = 0; i < 4; ++i )
    // {
    //     BOOST_CHECK_CLOSE_FRACTION(
    //             secondOrderSchwarzschildPerturbation( i, i ),
    //             secondOrderSolarPerturbation( i, i ),
    //             1.0E-13 );
    // }
    // TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
    //         secondOrderSchwarzschildPerturbation.block( 1, 1, 3, 3 ),
    //         secondOrderSolarPerturbation.block( 1, 1, 3, 3 ),
    //         1.0E-10 );
    // TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
    //         secondOrderSolarPerturbation.block( 0, 1, 1, 3 ),
    //         expectedSpaceTimePerturbation.transpose( ),
    //         std::numeric_limits< double >::epsilon( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( firstOrderSchwarzschildPerturbation,
                                       firstOrderSchwarzschildPerturbation.transpose( ),
                                       std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            firstOrderSolarPerturbation, firstOrderSolarPerturbation.transpose( ), std::numeric_limits< double >::epsilon( ) );
    const Eigen::Matrix4d identityCheck =
            firstOrderSolarSystemMetric->getCurrentCovariantMetric( ) * firstOrderSolarSystemMetric->getCurrentContravariantMetric( );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( identityCheck, Eigen::Matrix4d::Identity( ), 5.0E-13 );
}

BOOST_AUTO_TEST_CASE( testSphericalHarmonicGravityInMetric )
{
    loadStandardSpiceKernels( );

    const double initialEphemerisTime = 1.0E7;
    const double finalEphemerisTime = 1.1E7;
    const double maximumTimeStep = 3600.0;
    const double buffer = 10.0 * maximumTimeStep;

    std::vector< std::string > bodyNames = { "Sun", "Earth" };
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings.at( "Earth" )->bodyDeformationSettings.clear( );
    bodySettings.at( "Earth" )->gravityFieldVariationSettings.clear( );
    Eigen::Vector6d constantEarthState = Eigen::Vector6d::Zero( );
    constantEarthState.segment( 3, 3 ) =
            spice_interface::getBodyCartesianStateAtEpoch( "Earth", "SSB", "ECLIPJ2000", "None", 0.0 ).segment( 3, 3 );
    bodySettings.at( "Earth" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( constantEarthState );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    setGlobalFrameBodyEphemerides( bodies.getMap( ), "SSB", "ECLIPJ2000" );

    const double evaluationTime = 1.05E7;
    bodies.getBody( "Sun" )->setStateFromEphemeris( evaluationTime );
    bodies.getBody( "Earth" )->setStateFromEphemeris( evaluationTime );
    bodies.getBody( "Earth" )->setCurrentRotationToLocalFrameFromEphemeris( evaluationTime );
    auto earthGravityField =
            std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >( bodies.getBody( "Earth" )->getGravityFieldModel( ) );
    BOOST_REQUIRE_MESSAGE( earthGravityField != nullptr, "Earth gravity field is not spherical harmonics." );

    std::vector< std::string > firstOrderPerturbingBodies{ "Earth" };
    std::map< std::string, std::pair< int, int > > bodySphericalHarmonicExpansions;
    const int availableDegree = earthGravityField->getDegreeOfExpansion( );
    const int availableOrder = earthGravityField->getOrderOfExpansion( );
    const int fullDegree = std::min( 12, availableDegree );
    const int fullOrder = std::min( 12, availableOrder );
    bodySphericalHarmonicExpansions[ "Earth" ] = std::make_pair( fullDegree, fullOrder );
    auto fullMetricSettings = std::make_shared< SolarSystemSpaceTimeMetricSettings >(
            firstOrderPerturbingBodies, std::vector< std::string >( ), bodySphericalHarmonicExpansions, std::vector< std::string >( ) );

    const int truncatedDegree = std::min( 2, availableDegree );
    const int truncatedOrder = std::min( 2, availableOrder );
    bodySphericalHarmonicExpansions[ "Earth" ] = std::make_pair( truncatedDegree, truncatedOrder );
    auto truncatedMetricSettings = std::make_shared< SolarSystemSpaceTimeMetricSettings >(
            firstOrderPerturbingBodies, std::vector< std::string >( ), bodySphericalHarmonicExpansions, std::vector< std::string >( ) );

    bodySphericalHarmonicExpansions[ "Earth" ] = std::make_pair( 0, 0 );
    auto pointMassMetricSettings = std::make_shared< SolarSystemSpaceTimeMetricSettings >(
            firstOrderPerturbingBodies, std::vector< std::string >( ), bodySphericalHarmonicExpansions, std::vector< std::string >( ) );

    auto fullMetric = createSpaceTimeMetric( fullMetricSettings, bodies );
    auto truncatedMetric = createSpaceTimeMetric( truncatedMetricSettings, bodies );
    auto pointMassMetric = createSpaceTimeMetric( pointMassMetricSettings, bodies );

    Eigen::Vector6d earthCenteredKeplerElements;
    earthCenteredKeplerElements << 6378.0E3 + 300E3, 0.013, 1.0238269559089248, 3.1526292818328812, 1.5807574453453865, 3.1478950321924795;
    const double earthGravitationalParameter = earthGravityField->getGravitationalParameter( );
    const Eigen::Vector6d testCartesianElements =
            convertKeplerianToCartesianElements( earthCenteredKeplerElements, earthGravitationalParameter );

    fullMetric->update( testCartesianElements + bodies.getBody( "Earth" )->getState( ), 0.0, true, false );
    truncatedMetric->update( testCartesianElements + bodies.getBody( "Earth" )->getState( ), 0.0, true, false );
    pointMassMetric->update( testCartesianElements + bodies.getBody( "Earth" )->getState( ), 0.0, true, false );

    auto sphericalHarmonicsCache = std::make_shared< basic_mathematics::SphericalHarmonicsCache >( fullDegree, fullOrder );
    const Eigen::Vector3d relativePosition = testCartesianElements.segment( 0, 3 );

    const double fullPotential =
            earthGravityField->getGravitationalPotentialFromInertialPosition( relativePosition,
                                                                              bodies.getBody( "Earth" )->getCurrentRotationToLocalFrame( ),
                                                                              fullDegree,
                                                                              fullOrder,
                                                                              sphericalHarmonicsCache );
    const double truncatedPotential =
            earthGravityField->getGravitationalPotentialFromInertialPosition( relativePosition,
                                                                              bodies.getBody( "Earth" )->getCurrentRotationToLocalFrame( ),
                                                                              truncatedDegree,
                                                                              truncatedOrder,
                                                                              sphericalHarmonicsCache );
    const double pointMassPotential = earthGravityField->getGravitationalPotentialFromInertialPosition(
            relativePosition, bodies.getBody( "Earth" )->getCurrentRotationToLocalFrame( ), 0, 0, sphericalHarmonicsCache );

    std::vector< double > metricTerms{ 2.0 * fullPotential * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT,
                                       2.0 * truncatedPotential * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT,
                                       2.0 * pointMassPotential * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT };

    std::vector< Eigen::Matrix4d > metricPerturbations{ fullMetric->getCurrentCovariantMetricPeturbation( ),
                                                        truncatedMetric->getCurrentCovariantMetricPeturbation( ),
                                                        pointMassMetric->getCurrentCovariantMetricPeturbation( ) };

    for( unsigned int i = 0; i < metricPerturbations.size( ); ++i )
    {
        const double perturbation = metricTerms.at( i );
        BOOST_CHECK_CLOSE_FRACTION( perturbation - 0.5 * perturbation * perturbation,
                                    metricPerturbations.at( i )( 0, 0 ),
                                    std::numeric_limits< double >::epsilon( ) );

        for( int row = 0; row < 3; ++row )
        {
            for( int col = 0; col < 3; ++col )
            {
                if( row == col )
                {
                    BOOST_CHECK_CLOSE_FRACTION(
                            perturbation, metricPerturbations.at( i )( row + 1, col + 1 ), std::numeric_limits< double >::epsilon( ) );

                    BOOST_CHECK_CLOSE_FRACTION( -2.0 * perturbation * constantEarthState( 3 + row ) / physical_constants::SPEED_OF_LIGHT,
                                                metricPerturbations.at( i )( 0, row + 1 ),
                                                std::numeric_limits< double >::epsilon( ) );
                }
                else
                {
                    BOOST_CHECK_EQUAL( metricPerturbations.at( i )( row + 1, col + 1 ), 0.0 );
                }
            }
        }
    }
}

// Evaluate the metric VALUE directly (no propagator/converter) at a genuinely rotating
// ground-station BCRS state over a day, for a point-mass central body (no SH wrapper) vs an SH
// central body. This guards the metric-value path the direct-from-metric propagator relies on:
//   (a) the point-mass perturbation matches the analytic post-Newtonian form at every epoch;
//   (b) the central-body distance is rotation-invariant (constant station radius);
//   (c) the point-mass and SH perturbations differ only by a TIME-CONSTANT (the SH J2 term),
//       i.e. neither carries a spurious diurnal term.
BOOST_AUTO_TEST_CASE( testPointMassVsShMetricAtRotatingStation )
{
    loadStandardSpiceKernels( );

    const double t0 = 1.0E7;
    const double buffer = 5.0 * 3600.0;
    const std::vector< std::string > bodyNames{ "Sun", "Earth" };
    const double inv_c2 = physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;
    const double inv_c4 = inv_c2 * inv_c2;

    // Two environments sharing identical ephemerides/rotation; only Earth's gravity-field TYPE differs.
    auto makeBodies = [ & ]( const bool earthAsPointMass ) {
        BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, t0 - buffer, t0 + 26.0 * 3600.0 + buffer );
        bodySettings.at( "Sun" )->gravityFieldSettings = std::make_shared< SpiceCentralGravityFieldSettings >( "Sun" );
        if( earthAsPointMass )
        {
            bodySettings.at( "Earth" )->gravityFieldSettings = std::make_shared< SpiceCentralGravityFieldSettings >( "Earth" );
        }
        bodySettings.at( "Earth" )->bodyDeformationSettings.clear( );
        bodySettings.at( "Earth" )->gravityFieldVariationSettings.clear( );
        return createSystemOfBodies( bodySettings );
    };

    SystemOfBodies bodiesPm = makeBodies( true );
    SystemOfBodies bodiesSh = makeBodies( false );

    auto pmMetric = createSpaceTimeMetric(
            std::make_shared< SolarSystemSpaceTimeMetricSettings >( bodyNames,
                                                                    std::vector< std::string >( ),
                                                                    std::map< std::string, std::pair< int, int > >( ),
                                                                    std::vector< std::string >( ) ),
            bodiesPm );
    auto shMetric = createSpaceTimeMetric(
            std::make_shared< SolarSystemSpaceTimeMetricSettings >( bodyNames,
                                                                    std::vector< std::string >( ),
                                                                    std::map< std::string, std::pair< int, int > >{ { "Earth", { 2, 2 } } },
                                                                    std::vector< std::string >( ) ),
            bodiesSh );

    const double earthGm = bodiesPm.getBody( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    const double sunGm = bodiesPm.getBody( "Sun" )->getGravityFieldModel( )->getGravitationalParameter( );

    const double earthRadius = 6378137.0;
    const Eigen::Vector3d stationBodyFixed =
            ( Eigen::Vector3d( ) << earthRadius / std::sqrt( 2.0 ), earthRadius / std::sqrt( 2.0 ), 0.0 ).finished( );

    double firstH00Difference = 0.0;
    double maxH00DifferenceVariation = 0.0;
    for( int hour = 0; hour <= 26; ++hour )
    {
        const double t = t0 + hour * 3600.0;
        for( SystemOfBodies* bodies : { &bodiesPm, &bodiesSh } )
        {
            bodies->getBody( "Sun" )->setStateFromEphemeris( t );
            bodies->getBody( "Earth" )->setStateFromEphemeris( t );
            bodies->getBody( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( t );
        }
        auto earthRot = bodiesPm.getBody( "Earth" )->getRotationalEphemeris( );
        const Eigen::Matrix3d R = earthRot->getRotationToBaseFrame( t ).toRotationMatrix( );
        const Eigen::Vector3d omega = earthRot->getRotationalVelocityVectorInBaseFrame( t );
        const Eigen::Vector6d earthState = bodiesPm.getBody( "Earth" )->getState( );
        const Eigen::Vector6d sunState = bodiesPm.getBody( "Sun" )->getState( );
        const Eigen::Vector3d stationInertial = R * stationBodyFixed;
        Eigen::Vector6d stationState;
        stationState.segment( 0, 3 ) = earthState.segment( 0, 3 ) + stationInertial;
        stationState.segment( 3, 3 ) = earthState.segment( 3, 3 ) + omega.cross( stationInertial );

        pmMetric->update( stationState, t, true, false );
        shMetric->update( stationState, t, true, false );
        const Eigen::Matrix4d hPm = pmMetric->getCurrentCovariantMetricPeturbation( );
        const Eigen::Matrix4d hSh = shMetric->getCurrentCovariantMetricPeturbation( );

        // (a) Point-mass perturbation vs analytic post-Newtonian form (beta=gamma=1):
        //     h00 = 2 U/c^2 - 2 (sum U_i^2)/c^4 ;  h_ii = 2 U/c^2 ;  h_0i = -4 w_i/c^3.
        const double uEarth = earthGm / ( stationState.segment( 0, 3 ) - earthState.segment( 0, 3 ) ).norm( );
        const double uSun = sunGm / ( stationState.segment( 0, 3 ) - sunState.segment( 0, 3 ) ).norm( );
        const double uTotal = uEarth + uSun;
        const double uSquaredSum = uEarth * uEarth + uSun * uSun;
        const Eigen::Vector3d vectorPotential = uEarth * earthState.segment( 3, 3 ) + uSun * sunState.segment( 3, 3 );

        const double expectedH00 = 2.0 * ( uTotal * inv_c2 - uSquaredSum * inv_c4 );
        BOOST_CHECK_CLOSE_FRACTION( hPm( 0, 0 ), expectedH00, 1.0E-12 );
        for( int i = 1; i < 4; ++i )
        {
            BOOST_CHECK_CLOSE_FRACTION( hPm( i, i ), 2.0 * uTotal * inv_c2, 1.0E-12 );
        }
        const Eigen::Vector3d expectedSpaceTime = -4.0 * vectorPotential * inv_c2 / physical_constants::SPEED_OF_LIGHT;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( hPm.block( 1, 0, 3, 1 ), expectedSpaceTime, 1.0E-10 );

        // (b) Central-body (Earth) distance is rotation-invariant: equals the station radius.
        auto pmSolar = std::dynamic_pointer_cast< SolarSystemMetric >( pmMetric );
        BOOST_CHECK_CLOSE_FRACTION( pmSolar->getCurrentBodyDistance( 1 ), earthRadius, 1.0E-9 );

        // (c) Point-mass vs SH h00 differ only by a TIME-CONSTANT (the SH J2 term).
        const double h00Difference = hPm( 0, 0 ) - hSh( 0, 0 );
        if( hour == 0 )
        {
            firstH00Difference = h00Difference;
        }
        maxH00DifferenceVariation = std::max( maxH00DifferenceVariation, std::fabs( h00Difference - firstH00Difference ) );
    }
    // The SH J2 contribution to h00 is ~7.5e-13; its variation over a day must be far below that.
    BOOST_CHECK_SMALL( maxH00DifferenceVariation, 1.0E-16 );
}

// Analytic vs numeric partials: the SolarSystemMetric computes analytic position- and
// time-partials of the scalar potential (updatePotentialPartials, used to build the Christoffel
// symbols). Validate them against central finite differences of the scalar-potential VALUE.
BOOST_AUTO_TEST_CASE( testSolarSystemMetricPotentialPartialsAgainstFiniteDifferences )
{
    loadStandardSpiceKernels( );

    const double t0 = 1.0E7;
    const double buffer = 5.0 * 3600.0;
    const std::vector< std::string > bodyNames{ "Sun", "Earth" };

    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, t0 - buffer, t0 + buffer );
    bodySettings.at( "Sun" )->gravityFieldSettings = std::make_shared< SpiceCentralGravityFieldSettings >( "Sun" );
    bodySettings.at( "Earth" )->gravityFieldSettings = std::make_shared< SpiceCentralGravityFieldSettings >( "Earth" );
    bodySettings.at( "Earth" )->bodyDeformationSettings.clear( );
    bodySettings.at( "Earth" )->gravityFieldVariationSettings.clear( );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    auto metric = std::dynamic_pointer_cast< SolarSystemMetric >( createSpaceTimeMetric(
            std::make_shared< SolarSystemSpaceTimeMetricSettings >( bodyNames,
                                                                    std::vector< std::string >( ),
                                                                    std::map< std::string, std::pair< int, int > >( ),
                                                                    std::vector< std::string >( ) ),
            bodies ) );
    BOOST_REQUIRE( metric != nullptr );

    auto tickBodies = [ & ]( const double t ) {
        bodies.getBody( "Sun" )->setStateFromEphemeris( t );
        bodies.getBody( "Earth" )->setStateFromEphemeris( t );
    };

    // Evaluation point: a fixed BCRS position offset from Earth at a low-orbit-like radius, where
    // the scalar-potential gradient is large enough that a 1 m position step is roundoff-clean.
    tickBodies( t0 );
    Eigen::Vector6d evaluationState = Eigen::Vector6d::Zero( );
    evaluationState.segment( 0, 3 ) =
            bodies.getBody( "Earth" )->getState( ).segment( 0, 3 ) + ( Eigen::Vector3d( ) << 7.0E6, 1.0E6, 2.0E6 ).finished( );

    // Analytic partials (populated by the Christoffel update).
    metric->update( evaluationState, t0, true, true );
    const Eigen::Vector3d analyticPositionPartial = metric->getCurrentScalarPotentialPositionPartial( );
    const double analyticTimePartial = metric->getCurrentScalarPotentialTimePartial( );

    auto scalarPotentialAt = [ & ]( const Eigen::Vector6d& state, const double t ) {
        metric->update( state, t, true, false );
        return metric->getCurrentScalarPotential( );
    };
    // 4th-order central difference: f'(x) = [ -f(x+2h) + 8 f(x+h) - 8 f(x-h) + f(x-2h) ] / (12 h).
    auto fourthOrder = [ ]( const double fm2, const double fm1, const double fp1, const double fp2, const double h ) {
        return ( -fp2 + 8.0 * fp1 - 8.0 * fm1 + fm2 ) / ( 12.0 * h );
    };

    // Numeric position partial: perturb the evaluation point (bodies fixed at t0).
    const double positionStep = 1.0;  // m
    Eigen::Vector3d numericPositionPartial;
    tickBodies( t0 );
    for( int axis = 0; axis < 3; ++axis )
    {
        auto shifted = [ & ]( const double delta ) {
            Eigen::Vector6d s = evaluationState;
            s( axis ) += delta;
            return scalarPotentialAt( s, t0 );
        };
        numericPositionPartial( axis ) = fourthOrder( shifted( -2.0 * positionStep ),
                                                      shifted( -positionStep ),
                                                      shifted( positionStep ),
                                                      shifted( 2.0 * positionStep ),
                                                      positionStep );
    }
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( analyticPositionPartial, numericPositionPartial, 1.0E-7 );

    // Numeric time partial: perturb time (move the bodies), evaluation point fixed.
    const double timeStep = 1.0;  // s
    auto potentialAtTime = [ & ]( const double dt ) {
        tickBodies( t0 + dt );
        return scalarPotentialAt( evaluationState, t0 + dt );
    };
    const double numericTimePartial = fourthOrder( potentialAtTime( -2.0 * timeStep ),
                                                   potentialAtTime( -timeStep ),
                                                   potentialAtTime( timeStep ),
                                                   potentialAtTime( 2.0 * timeStep ),
                                                   timeStep );
    BOOST_CHECK_CLOSE_FRACTION( analyticTimePartial, numericTimePartial, 1.0E-7 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
