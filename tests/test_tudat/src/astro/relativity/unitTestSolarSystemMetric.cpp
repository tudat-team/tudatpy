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

#include <limits>
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
    BodyListSettings bodySettings = getDefaultBodySettings(
            bodyNames,
            initialEphemerisTime - buffer,
            finalEphemerisTime + buffer );
    bodySettings.at( "Earth" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( Eigen::Vector6d::Zero( ) );
    bodySettings.at( "Earth" )->gravityFieldVariationSettings.clear( );
    auto bodies = createSystemOfBodies( bodySettings );
    std::cerr << "Bodies created" << std::endl;
    setGlobalFrameBodyEphemerides( bodies.getMap( ), "SSB", "ECLIPJ2000" );
    std::cerr << "Frame set" << std::endl;
    bodies.getBody( "Earth" )->setStateFromEphemeris( evaluationTime );
    std::cerr << "Earth initialized" << std::endl;

    auto firstOrderSchwarzschildSettings =
            std::make_shared< SchwardschildSpaceTimeMetricSettings >( "Earth", ppnParameterSet, false );
    [[maybe_unused]] auto secondOrderSchwarzschildSettings =
            std::make_shared< SchwardschildSpaceTimeMetricSettings >( "Earth", ppnParameterSet, true );

    auto firstOrderSolarSystemSettings =
            std::make_shared< SolarSystemSpaceTimeMetricSettings >(
                    std::vector< std::string >{ "Earth" },
                    std::vector< std::string >( ),
                    std::map< std::string, std::pair< int, int > >( ),
                    std::vector< std::string >( ),
                    ppnParameterSet );
    // Second-order solar-system terms currently unimplemented; re-enable once metric supports them.
    // auto secondOrderSolarSystemSettings =
    //         std::make_shared< SolarSystemSpaceTimeMetricSettings >(
    //                 std::vector< std::string >( ),
    //                 std::vector< std::string >{ "Earth" },
    //                 std::map< std::string, std::pair< int, int > >( ),
    //                 std::vector< std::string >( ),
    //                 ppnParameterSet );

    auto firstOrderSchwarzschildMetric = createSpaceTimeMetric( firstOrderSchwarzschildSettings, bodies );
    [[maybe_unused]] auto secondOrderSchwarzschildMetric = createSpaceTimeMetric( secondOrderSchwarzschildSettings, bodies );
    auto firstOrderSolarSystemMetric = createSpaceTimeMetric( firstOrderSolarSystemSettings, bodies );

    Eigen::Vector6d keplerianElements;
    keplerianElements << 6378.0E3 + 400.0E3, 0.013, 1.0238269559089248,
            3.1526292818328812, 1.5807574453453865, 3.1478950321924795;

    const double earthGravitationalParameter =
            bodies.getBody( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    const Eigen::Vector6d cartesianState = convertKeplerianToCartesianElements(
            keplerianElements, earthGravitationalParameter );

    firstOrderSchwarzschildMetric->update( cartesianState, 0.0, true, true );
    secondOrderSchwarzschildMetric->update( cartesianState, 0.0, true, true );

    const Eigen::Vector6d barycentricState = cartesianState + bodies.getBody( "Earth" )->getState( );
    firstOrderSolarSystemMetric->update( barycentricState, 0.0, true, true );

    const Eigen::Matrix4d firstOrderSchwarzschildPerturbation =
            firstOrderSchwarzschildMetric->getCurrentCovariantMetricPeturbation( );
    const Eigen::Matrix4d firstOrderSolarPerturbation =
            firstOrderSolarSystemMetric->getCurrentCovariantMetricPeturbation( );

    BOOST_CHECK_CLOSE_FRACTION(
            firstOrderSchwarzschildPerturbation( 0, 0 ),
            firstOrderSolarPerturbation( 0, 0 ),
            1.0E-13 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            firstOrderSchwarzschildPerturbation.block( 1, 1, 3, 3 ),
            firstOrderSolarPerturbation.block( 1, 1, 3, 3 ),
            1.0E-13 );

    const Eigen::Vector3d earthVelocity = bodies.getBody( "Earth" )->getState( ).segment( 3, 3 );
    const Eigen::Vector3d expectedSpaceTimePerturbation =
            -2.0 * earthVelocity * firstOrderSolarPerturbation( 1, 1 ) / physical_constants::SPEED_OF_LIGHT;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            firstOrderSolarPerturbation.block( 0, 1, 1, 3 ),
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

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            firstOrderSchwarzschildPerturbation,
            firstOrderSchwarzschildPerturbation.transpose( ),
            std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            firstOrderSolarPerturbation,
            firstOrderSolarPerturbation.transpose( ),
            std::numeric_limits< double >::epsilon( ) );
    const Eigen::Matrix4d identityCheck =
            firstOrderSolarSystemMetric->getCurrentCovariantMetric( ) *
            firstOrderSolarSystemMetric->getCurrentContravariantMetric( );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            identityCheck,
            Eigen::Matrix4d::Identity( ),
            5.0E-13 );
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
            spice_interface::getBodyCartesianStateAtEpoch( "Earth", "SSB" , "ECLIPJ2000" , "None", 0.0 ).segment( 3, 3 );
    bodySettings.at( "Earth" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( constantEarthState );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    setGlobalFrameBodyEphemerides( bodies.getMap( ), "SSB", "ECLIPJ2000" );

    const double evaluationTime = 1.05E7;
    bodies.getBody( "Sun" )->setStateFromEphemeris( evaluationTime );
    bodies.getBody( "Earth" )->setStateFromEphemeris( evaluationTime );
    bodies.getBody( "Earth" )->setCurrentRotationToLocalFrameFromEphemeris( evaluationTime );
    auto earthGravityField =
            std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                    bodies.getBody( "Earth" )->getGravityFieldModel( ) );
    BOOST_REQUIRE_MESSAGE( earthGravityField != nullptr, "Earth gravity field is not spherical harmonics." );

    std::vector< std::string > firstOrderPerturbingBodies{ "Earth" };
    std::map< std::string, std::pair< int, int > > bodySphericalHarmonicExpansions;
    const int availableDegree = earthGravityField->getDegreeOfExpansion( );
    const int availableOrder = earthGravityField->getOrderOfExpansion( );
    const int fullDegree = std::min( 12, availableDegree );
    const int fullOrder = std::min( 12, availableOrder );
    auto ppnSet = std::make_shared< relativity::PPNParameterSet >( 1.0, 1.0 );

    bodySphericalHarmonicExpansions[ "Earth" ] = std::make_pair( fullDegree, fullOrder );
    auto fullMetricSettings = std::make_shared< SolarSystemSpaceTimeMetricSettings >(
            firstOrderPerturbingBodies, std::vector< std::string >( ), bodySphericalHarmonicExpansions, std::vector< std::string >( ), ppnSet );

    const int truncatedDegree = std::min( 2, availableDegree );
    const int truncatedOrder = std::min( 2, availableOrder );
    bodySphericalHarmonicExpansions[ "Earth" ] = std::make_pair( truncatedDegree, truncatedOrder );
    auto truncatedMetricSettings = std::make_shared< SolarSystemSpaceTimeMetricSettings >(
            firstOrderPerturbingBodies, std::vector< std::string >( ), bodySphericalHarmonicExpansions, std::vector< std::string >( ), ppnSet );

    bodySphericalHarmonicExpansions[ "Earth" ] = std::make_pair( 0, 0 );
    auto pointMassMetricSettings = std::make_shared< SolarSystemSpaceTimeMetricSettings >(
            firstOrderPerturbingBodies, std::vector< std::string >( ), bodySphericalHarmonicExpansions, std::vector< std::string >( ), ppnSet );

    auto fullMetric = createSpaceTimeMetric( fullMetricSettings, bodies );
    auto truncatedMetric = createSpaceTimeMetric( truncatedMetricSettings, bodies );
    auto pointMassMetric = createSpaceTimeMetric( pointMassMetricSettings, bodies );

    Eigen::Vector6d earthCenteredKeplerElements;
    earthCenteredKeplerElements << 6378.0E3 + 300E3, 0.013, 1.0238269559089248,
            3.1526292818328812, 1.5807574453453865, 3.1478950321924795;
    const double earthGravitationalParameter = earthGravityField->getGravitationalParameter( );
    const Eigen::Vector6d testCartesianElements = convertKeplerianToCartesianElements(
            earthCenteredKeplerElements, earthGravitationalParameter );

    fullMetric->update( testCartesianElements + bodies.getBody( "Earth" )->getState( ), 0.0, true, false );
    truncatedMetric->update( testCartesianElements + bodies.getBody( "Earth" )->getState( ), 0.0, true, false );
    pointMassMetric->update( testCartesianElements + bodies.getBody( "Earth" )->getState( ), 0.0, true, false );

    auto legendreCache = std::make_unique< basic_mathematics::LegendreCache >( fullDegree, fullOrder );
    const Eigen::Vector3d relativePosition = testCartesianElements.segment( 0, 3 );

    const double fullPotential = earthGravityField->getGravitationalPotentialFromInertialPosition(
            relativePosition, bodies.getBody( "Earth" )->getCurrentRotationToLocalFrame( ), fullDegree, fullOrder, legendreCache.get( ) );
    const double truncatedPotential = earthGravityField->getGravitationalPotentialFromInertialPosition(
            relativePosition, bodies.getBody( "Earth" )->getCurrentRotationToLocalFrame( ), truncatedDegree, truncatedOrder, legendreCache.get( ) );
    const double pointMassPotential = earthGravityField->getGravitationalPotentialFromInertialPosition(
            relativePosition, bodies.getBody( "Earth" )->getCurrentRotationToLocalFrame( ), 0, 0, legendreCache.get( ) );

    std::vector< double > metricTerms{
        2.0 * fullPotential * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT,
        2.0 * truncatedPotential * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT,
        2.0 * pointMassPotential * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT };

    std::vector< Eigen::Matrix4d > metricPerturbations{
        fullMetric->getCurrentCovariantMetricPeturbation( ),
        truncatedMetric->getCurrentCovariantMetricPeturbation( ),
        pointMassMetric->getCurrentCovariantMetricPeturbation( ) };

    for ( unsigned int i = 0; i < metricPerturbations.size( ); ++i )
    {
        const double perturbation = metricTerms.at( i );
        BOOST_CHECK_CLOSE_FRACTION(
                perturbation - 0.5 * perturbation * perturbation,
                metricPerturbations.at( i )( 0, 0 ),
                std::numeric_limits< double >::epsilon( ) );

        for ( int row = 0; row < 3; ++row )
        {
            for ( int col = 0; col < 3; ++col )
            {
                if ( row == col )
                {
                    BOOST_CHECK_CLOSE_FRACTION(
                            perturbation,
                            metricPerturbations.at( i )( row + 1, col + 1 ),
                            std::numeric_limits< double >::epsilon( ) );

                    BOOST_CHECK_CLOSE_FRACTION(
                            -2.0 * perturbation * constantEarthState( 3 + row ) / physical_constants::SPEED_OF_LIGHT,
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

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
