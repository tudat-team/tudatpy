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

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/simulation/simulation.h"
#include "tudat/math/basic/leastSquaresEstimation.h"
#include "tudat/io/basicInputOutput.h"

namespace tudat
{

namespace unit_tests
{

using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::gravitation;
using namespace tudat::basic_astrodynamics;
using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::numerical_integrators;


BOOST_AUTO_TEST_SUITE( test_mutual_spherical_harmonic_gravity )

double computeSemiMajorAxisRateDueToTideRaisedOnPlanet(
        const double k2QRatio, const double satelliteWrtPlanetMassRatio, const double referenceRadiusWrtSemiMajorAxisRatio,
        const double semiMajorAxis, const double meanMotion )
{
    return  3.0 * k2QRatio * satelliteWrtPlanetMassRatio * std::pow( referenceRadiusWrtSemiMajorAxisRatio, 5.0 ) *
            semiMajorAxis * meanMotion;
}

double computeEccentricityRateDueToTideRaisedOnPlanet(
        const double k2QRatio, const double satelliteWrtPlanetMassRatio, const double referenceRadiusWrtSemiMajorAxisRatio,
        const double eccentricity, const double meanMotion )
{
    return 57.0 / 8.0 * k2QRatio * satelliteWrtPlanetMassRatio * std::pow( referenceRadiusWrtSemiMajorAxisRatio, 5.0 ) *
            eccentricity * meanMotion;
}


std::pair< double, double > computeKeplerElementRatesDueToDissipation(
        const SystemOfBodies& bodies, const std::string& satelliteToPropagate, const bool usePlanetDissipation,
        const double k2LoveNumber, const double tidalTimeLag, const double initialTime, const double finalTime,
        double& meanMotion, Eigen::Vector6d& intialKeplerElements  )
{
    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    bodiesToPropagate.push_back( satelliteToPropagate );
    centralBodies.push_back( "Jupiter" );

    // Define propagation settings.
    const double fixedStepSize = 450.0;
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< RungeKuttaVariableStepSizeSettings< > >
            ( 0.0, fixedStepSize,
              rungeKuttaFehlberg78, fixedStepSize, fixedStepSize, 1.0, 1.0);

    std::map< double, Eigen::VectorXd > integrationResultWithDissipation;
    std::map< double, Eigen::VectorXd > integrationResultWithDissipationKepler;

    std::map< double, double > semiMajorAxes, eccentricities;
    {
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfIo;
        accelerationsOfIo[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >(
                                                      basic_astrodynamics::point_mass_gravity ) );
        accelerationsOfIo[ "Jupiter" ].push_back( std::make_shared< DirectTidalDissipationAccelerationSettings >(
                                                      k2LoveNumber, tidalTimeLag, false, usePlanetDissipation ) );
        accelerationMap[ satelliteToPropagate ] = accelerationsOfIo;
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodies, accelerationMap, bodiesToPropagate, centralBodies );


        // Save dependent variables
        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesToSave;
        if( usePlanetDissipation )
        {
            dependentVariablesToSave.push_back(
                        std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                            direct_tidal_dissipation_in_central_body_acceleration, satelliteToPropagate, "Jupiter" ) );
        }
        else
        {
            dependentVariablesToSave.push_back(
                        std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                            direct_tidal_dissipation_in_orbiting_body_acceleration, satelliteToPropagate, "Jupiter" ) );
        }

        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, getInitialStatesOfBodies(
                      bodiesToPropagate, centralBodies, bodies, initialTime ), initialTime, integratorSettings,
                  std::make_shared< PropagationTimeTerminationSettings >( finalTime ), cowell,
                  dependentVariablesToSave );


        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodies, propagatorSettings );
        integrationResultWithDissipation = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

        for( std::map< double, Eigen::VectorXd >::const_iterator mapIterator = integrationResultWithDissipation.begin( );
             mapIterator != integrationResultWithDissipation.end( ); mapIterator++ )
        {
            integrationResultWithDissipationKepler[ mapIterator->first ] =
                    orbital_element_conversions::convertCartesianToKeplerianElements< double >(
                        mapIterator->second, getBodyGravitationalParameter( "Jupiter" ) +
                        getBodyGravitationalParameter( satelliteToPropagate ) );
            semiMajorAxes[ mapIterator->first ] = integrationResultWithDissipationKepler[ mapIterator->first ]( 0 ) -
                    integrationResultWithDissipationKepler.begin( )->second( 0 );
            eccentricities[ mapIterator->first ] = integrationResultWithDissipationKepler[ mapIterator->first ]( 1 ) -
                    integrationResultWithDissipationKepler.begin( )->second( 1 );

        }
    }

    //    input_output::writeDataMapToTextFile( integrationResultWithDissipationKepler,
    //                                          "keplerElements_"  + std::to_string( usePlanetDissipation ) +
    //                                          satelliteToPropagate + ".dat" );

    std::vector< double > semiMajorAxisFit = linear_algebra::getLeastSquaresPolynomialFit(
                semiMajorAxes, { 0, 1 } );
    std::vector< double > eccentricityFit = linear_algebra::getLeastSquaresPolynomialFit(
                eccentricities, { 0, 1 } );

    intialKeplerElements = integrationResultWithDissipationKepler.begin( )->second;
    meanMotion = basic_astrodynamics::computeKeplerMeanMotion(
                intialKeplerElements( 0 ), getBodyGravitationalParameter( "Jupiter" ) +
                getBodyGravitationalParameter( satelliteToPropagate ) );

    return std::make_pair( semiMajorAxisFit.at( 1 ), eccentricityFit.at( 1 ) );
}

BOOST_AUTO_TEST_CASE( testTidalDissipationInPlanetAndSatellite )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create list of bodies to create.
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Io" );
    bodyNames.push_back( "Europa" );
    bodyNames.push_back( "Ganymede" );
    //bodyNames.push_back( "Callisto" );
    bodyNames.push_back( "Sun" );

    // Specify initial time
    double initialTime = 0.0;
    double finalTime = 1.0 * physical_constants::JULIAN_YEAR;

    // Get body settings.
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames, initialTime - 86400.0, finalTime + 86400.0 );

    std::vector< std::string > galileanSatellites = { "Io", "Europa" };

    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
    cosineCoefficients( 0, 0 ) = 1.0;
    Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
    bodySettings.at( "Jupiter" )->gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >
            ( getBodyGravitationalParameter( "Jupiter" ), getAverageRadius( "Jupiter" ),
              cosineCoefficients, sineCoefficients, "IAU_Jupiter" );
    bodySettings.at( "Jupiter" )->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                "ECLIPJ2000", "IAU_Jupiter", Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ),
                0.0, 2.0 * mathematical_constants::PI / ( 9.925 * 3600.0 ) );

    for( unsigned int i = 0; i < galileanSatellites.size( ); i++ )
    {
        bodySettings.at( galileanSatellites.at( i ) )->gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >
                ( getBodyGravitationalParameter( galileanSatellites.at( i ) ), getAverageRadius( galileanSatellites.at( i ) ),
                  cosineCoefficients, sineCoefficients, "IAU_" + galileanSatellites.at( i )  );
    }

    bodySettings.at( "Io" )->ephemerisSettings = std::make_shared< KeplerEphemerisSettings >(
                ( Eigen::Vector6d( ) << 1.0 * 421.8E6, 1.0 * 0.004, 0.0, 0.0, 0.0, 0.0 ).finished( ), 0.0,
                getBodyGravitationalParameter( "Jupiter" ) + getBodyGravitationalParameter( "Io" ), "Jupiter", "ECLIPJ2000" );
    bodySettings.at( "Europa" )->ephemerisSettings = std::make_shared< KeplerEphemerisSettings >(
                ( Eigen::Vector6d( ) << 671.1E6, 0.009, 0.0, 0.0, 0.0, 0.0 ).finished( ), 0.0,
                getBodyGravitationalParameter( "Jupiter" ) + getBodyGravitationalParameter( "Europa" ), "Jupiter", "ECLIPJ2000" );
    bodySettings.at( "Ganymede" )->ephemerisSettings = std::make_shared< KeplerEphemerisSettings >(
                ( Eigen::Vector6d( ) << 1070.400E6, 0.0013, 0.0, 0.0, 0.0, 0.0 ).finished( ), 0.0,
                getBodyGravitationalParameter( "Jupiter" ) + getBodyGravitationalParameter( "Ganymede" ), "Jupiter", "ECLIPJ2000" );
    //    bodySettings[ "Callisto" ]->ephemerisSettings = std::make_shared< KeplerEphemerisSettings >(
    //                ( Eigen::Vector6d( ) << 1882.700E6, 0.0074, 0.0, 0.0, 0.0, 0.0 ).finished( ), 0.0,
    //                getBodyGravitationalParameter( "Jupiter" ) + getBodyGravitationalParameter( "Callisto" ), "Jupiter", "ECLIPJ2000" );

    // Create bodies needed in simulation
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    

    // Define propagation settings.
    double jupiterLoveNumber = 0.1;
    double jupiterTimeLag = 100.0;

    for( unsigned int i = 0; i < galileanSatellites.size( ); i++ )
    {
        Eigen::Vector6d intialKeplerElements;
        double meanMotion;
        std::pair< double, double > elementRates = computeKeplerElementRatesDueToDissipation(
                    bodies, galileanSatellites.at( i ), true, jupiterLoveNumber, jupiterTimeLag, initialTime, finalTime, meanMotion, intialKeplerElements );

        double orbitalPeriod = 2.0 * mathematical_constants::PI / meanMotion;
        double jupiterForcingTime = 9.925 * 3600.0 * orbitalPeriod / ( 2.0 * std::fabs( orbitalPeriod - 9.925 * 3600.0 ) );
        double jupiterQualityFactor = 1.0 / std::sin( 2.0 * mathematical_constants::PI * jupiterTimeLag / ( jupiterForcingTime ) );

        double theoreticalSemiMajorAxisRateFromJupiterTide = computeSemiMajorAxisRateDueToTideRaisedOnPlanet(
                    jupiterLoveNumber / jupiterQualityFactor, getBodyGravitationalParameter( galileanSatellites.at( i ) ) /
                    getBodyGravitationalParameter( "Jupiter" ), getAverageRadius( "Jupiter" ) / intialKeplerElements( 0 ),
                    intialKeplerElements( 0 ), meanMotion );
        double theoreticaEccentricityRateFromJupiterTide = computeEccentricityRateDueToTideRaisedOnPlanet(
                    jupiterLoveNumber / jupiterQualityFactor, getBodyGravitationalParameter( galileanSatellites.at( i ) ) /
                    getBodyGravitationalParameter( "Jupiter" ), getAverageRadius( "Jupiter" ) / intialKeplerElements( 0 ),
                    intialKeplerElements( 1 ), meanMotion );

        //BOOST_CHECK_CLOSE_FRACTION( elementRates.first, theoreticalSemiMajorAxisRateFromJupiterTide, 2.0E-3 );
        BOOST_CHECK( std::fabs(elementRates.first - theoreticalSemiMajorAxisRateFromJupiterTide) < std::fabs(std::min(elementRates.first, theoreticalSemiMajorAxisRateFromJupiterTide) * 2.0E-3 ));
        //BOOST_CHECK_CLOSE_FRACTION( elementRates.second, theoreticaEccentricityRateFromJupiterTide, 1.0E-1 );
        BOOST_CHECK( std::fabs(elementRates.second - theoreticaEccentricityRateFromJupiterTide) < std::fabs(std::min(elementRates.second, theoreticaEccentricityRateFromJupiterTide) * 1.0E-1 ));
    }


    std::map< double, Eigen::VectorXd > integrationResultWithDissipationInIo;
    std::map< double, Eigen::VectorXd > integrationResultWithDissipationInIoKepler;

    double satelliteLoveNumber = 1.0;
    double satelliteTimeLag = 1000.0;

    for( unsigned int i = 0; i < galileanSatellites.size( ); i++ )
    {
        Eigen::Vector6d intialKeplerElements;
        double meanMotion;
        std::pair< double, double > elementRates = computeKeplerElementRatesDueToDissipation(
                    bodies, galileanSatellites.at( i ), false, satelliteLoveNumber, satelliteTimeLag, initialTime, finalTime, meanMotion, intialKeplerElements );

        double orbitalPeriod = 2.0 * mathematical_constants::PI / meanMotion;
        double satelliteQualityFactor = 1.0 / std::sin( 2.0 * mathematical_constants::PI * satelliteTimeLag / orbitalPeriod );

        double theoreticalSemiMajorAxisRateFromIoTide =  -21.0 * satelliteLoveNumber / satelliteQualityFactor * getBodyGravitationalParameter( "Jupiter" ) /
                getBodyGravitationalParameter( galileanSatellites.at( i ) ) * std::pow( getAverageRadius( galileanSatellites.at( i ) ) / intialKeplerElements( 0 ), 5.0 ) *
                intialKeplerElements( 0 ) * meanMotion * intialKeplerElements( 1 ) * intialKeplerElements( 1 );
        double theoreticaEccentricityRateFromIoTide =  - 21.0 / 2.0 * satelliteLoveNumber / satelliteQualityFactor * getBodyGravitationalParameter( "Jupiter" ) /
                getBodyGravitationalParameter( galileanSatellites.at( i ) ) * std::pow( getAverageRadius( galileanSatellites.at( i ) ) / intialKeplerElements( 0 ), 5.0 ) *
                intialKeplerElements( 1 ) * meanMotion;

        // Increase tolerance for more distance moons
        double toleranceMultiplier = 1.0;
        if( i == 1 )
        {
            toleranceMultiplier = 5.0;
        }
        else if( i == 2 )
        {
            toleranceMultiplier = 20.0;
        }

        //BOOST_CHECK_CLOSE_FRACTION( elementRates.first, theoreticalSemiMajorAxisRateFromIoTide, toleranceMultiplier * 1.0E-3 );
        BOOST_CHECK( std::fabs(elementRates.first - theoreticalSemiMajorAxisRateFromIoTide) < std::fabs(std::min(elementRates.first, theoreticalSemiMajorAxisRateFromIoTide) * toleranceMultiplier * 1.0E-3 ));
        //BOOST_CHECK_CLOSE_FRACTION( elementRates.second, theoreticaEccentricityRateFromIoTide, toleranceMultiplier * 1.0E-3 );
        BOOST_CHECK( std::fabs(elementRates.second - theoreticaEccentricityRateFromIoTide) < std::fabs(std::min(elementRates.second, theoreticaEccentricityRateFromIoTide) * toleranceMultiplier * 1.0E-3 ));

        // Artificially increase time lag to make effect observable over integration tiem of 1 year.
        satelliteTimeLag *= 5.0;
    }



}

BOOST_AUTO_TEST_SUITE_END( )

}

}
