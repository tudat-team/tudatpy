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

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Mathematics/BasicMathematics/leastSquaresEstimation.h"

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

BOOST_AUTO_TEST_CASE( testMutualSphericalHarmonicGravity )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create list of bodies to create.
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Io" );
    bodyNames.push_back( "Europa" );
    bodyNames.push_back( "Sun" );

    // Specify initial time
    double initialTime = 0.0;
    double finalTime = 10.0 * physical_constants::JULIAN_YEAR;

    // Get body settings.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames, initialTime - 86400.0, finalTime + 86400.0 );


    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
    cosineCoefficients( 0, 0 ) = 1.0;
    Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
    bodySettings[ "Jupiter" ]->gravityFieldSettings = boost::make_shared< SphericalHarmonicsGravityFieldSettings >
            ( getBodyGravitationalParameter( "Jupiter" ), getAverageRadius( "Jupiter" ),
              cosineCoefficients, sineCoefficients, "IAU_Jupiter" );
    bodySettings[ "Jupiter" ]->rotationModelSettings = boost::make_shared< SimpleRotationModelSettings >(
                "ECLIPJ2000", "IAU_Jupiter", Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ),
                0.0, 2.0 * mathematical_constants::PI / ( 9.925 * 3600.0 ) );

    bodySettings[ "Io" ]->gravityFieldSettings = boost::make_shared< SphericalHarmonicsGravityFieldSettings >
            ( getBodyGravitationalParameter( "Io" ), getAverageRadius( "Io" ),
              cosineCoefficients, sineCoefficients, "IAU_Io" );

    bodySettings[ "Io" ]->ephemerisSettings = boost::make_shared< KeplerEphemerisSettings >(
                ( Eigen::Vector6d( )<< 421.8E6, 0.004, 0.0, 0.0, 0.0, 0.0 ).finished( ), 0.0,
                getBodyGravitationalParameter( "Jupiter" ) + getBodyGravitationalParameter( "Io" ), "Jupiter", "ECLIPJ2000" );
    bodySettings[ "Europa" ]->ephemerisSettings = boost::make_shared< KeplerEphemerisSettings >(
                ( Eigen::Vector6d( )<< 671.1E6, 0.009, 0.0, 0.0, 0.0, 0.0 ).finished( ), 0.0,
                getBodyGravitationalParameter( "Jupiter" ) + getBodyGravitationalParameter( "Europa" ), "Jupiter", "ECLIPJ2000" );

    // Create bodies needed in simulation
    NamedBodyMap bodyMap = createBodies( bodySettings );
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    bodiesToPropagate.push_back( "Io" );
    centralBodies.push_back( "Jupiter" );

    // Define propagation settings.
    double jupiterLoveNumber = 0.1;
    double jupiterTimeLag = 100.0;

    const double fixedStepSize = 1800.0;
    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
            ( rungeKuttaVariableStepSize, 0.0, fixedStepSize,
              RungeKuttaCoefficients::rungeKuttaFehlberg78, 900.0, 900.0, 1.0, 1.0);

    //    std::map< double, Eigen::VectorXd > integrationResultWithDissipation;
    //    std::map< double, Eigen::VectorXd > integrationResultWithDissipationKepler;

    //    std::map< double, double > semiMajorAxes, eccentricities;
    //    {
    //        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfIo;
    //        accelerationsOfIo[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >(
    //                                                      basic_astrodynamics::central_gravity ) );
    //        accelerationsOfIo[ "Jupiter" ].push_back( boost::make_shared< DirectTidalDissipationAccelerationSettings >(
    //                                                      jupiterLoveNumber, jupiterTimeLag, false ) );
    //        accelerationMap[ "Io" ] = accelerationsOfIo;
    //        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
    //                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    //        boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
    //                boost::make_shared< TranslationalStatePropagatorSettings< double > >
    //                ( centralBodies, accelerationModelMap, bodiesToPropagate, getInitialStatesOfBodies(
    //                      bodiesToPropagate, centralBodies, bodyMap, initialTime ), finalTime );


    //        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    //        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //        // Create simulation object and propagate dynamics.
    //        SingleArcDynamicsSimulator< > dynamicsSimulator(
    //                    bodyMap, integratorSettings, propagatorSettings, true, false, false );
    //        integrationResultWithDissipation = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    //        for( std::map< double, Eigen::VectorXd >::const_iterator mapIterator = integrationResultWithDissipation.begin( );
    //             mapIterator != integrationResultWithDissipation.end( ); mapIterator++ )
    //        {
    //            integrationResultWithDissipationKepler[ mapIterator->first ] =
    //                    orbital_element_conversions::convertCartesianToKeplerianElements< double >(
    //                        mapIterator->second, getBodyGravitationalParameter( "Jupiter" ) +
    //                        getBodyGravitationalParameter( "Io" ) );
    //            semiMajorAxes[ mapIterator->first ] = integrationResultWithDissipationKepler[ mapIterator->first ]( 0 ) -
    //                    integrationResultWithDissipationKepler.begin( )->second( 0 );
    //            eccentricities[ mapIterator->first ] = integrationResultWithDissipationKepler[ mapIterator->first ]( 1 ) -
    //                    integrationResultWithDissipationKepler.begin( )->second( 1 );

    //        }
    //    }

    //    std::vector< double > semiMajorAxisFit = linear_algebra::getLeastSquaresPolynomialFit(
    //                semiMajorAxes, boost::assign::list_of( 0 )( 1 ) );
    //    std::vector< double > eccentricityFit = linear_algebra::getLeastSquaresPolynomialFit(
    //                eccentricities, boost::assign::list_of( 0 )( 1 ) );

    //    Eigen::Vector6d intialKeplerElements = integrationResultWithDissipationKepler.begin( )->second;
    //    double meanMotion = basic_astrodynamics::computeKeplerMeanMotion(
    //                intialKeplerElements( 0 ), getBodyGravitationalParameter( "Jupiter" ) + getBodyGravitationalParameter( "Io" ) );
    //    double orbitalPeriod = 2.0 * mathematical_constants::PI / meanMotion;
    //    double jupiterForcingTime = 9.925 * 3600.0 * orbitalPeriod / ( 2.0 * std::fabs( orbitalPeriod - 9.925 * 3600.0 ) );
    //    double jupiterQualityFactor = 1.0 / std::sin( 2.0 * mathematical_constants::PI * jupiterTimeLag / ( jupiterForcingTime ) );

    //    double theoreticalSemiMajorAxisRateFromJupiterTide =  3.0 * jupiterLoveNumber / jupiterQualityFactor * getBodyGravitationalParameter( "Io" ) /
    //            getBodyGravitationalParameter( "Jupiter" ) * std::pow( getAverageRadius( "Jupiter" ) / intialKeplerElements( 0 ), 5.0 ) *
    //            intialKeplerElements( 0 ) * meanMotion;
    //    double theoreticaEccentricityRateFromJupiterTide =  57.0 / 8.0 * jupiterLoveNumber / jupiterQualityFactor * getBodyGravitationalParameter( "Io" ) /
    //            getBodyGravitationalParameter( "Jupiter" ) * std::pow( getAverageRadius( "Jupiter" ) / intialKeplerElements( 0 ), 5.0 ) *
    //            intialKeplerElements( 1 ) * meanMotion;

    //    BOOST_CHECK_CLOSE_FRACTION( semiMajorAxisFit.at( 1 ), theoreticalSemiMajorAxisRateFromJupiterTide, 2.0E-3 );
    //    BOOST_CHECK_CLOSE_FRACTION( eccentricityFit.at( 1 ), theoreticaEccentricityRateFromJupiterTide, 1.0E-1 );


    std::map< double, Eigen::VectorXd > integrationResultWithDissipationInIo;
    std::map< double, Eigen::VectorXd > integrationResultWithDissipationInIoKepler;

    double ioLoveNumber = 1.0;
    double ioTimeLag = 1000.0;

    std::map< double, double > semiMajorAxesIoTide, eccentricitiesIoTide;
    {
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfIo;
        accelerationsOfIo[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >(
                                                      basic_astrodynamics::central_gravity ) );
        accelerationsOfIo[ "Jupiter" ].push_back( boost::make_shared< DirectTidalDissipationAccelerationSettings >(
                                                      ioLoveNumber, ioTimeLag, false, false ) );
        accelerationMap[ "Io" ] = accelerationsOfIo;
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

        boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                boost::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, getInitialStatesOfBodies(
                      bodiesToPropagate, centralBodies, bodyMap, initialTime ), finalTime );


        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodyMap, integratorSettings, propagatorSettings, true, false, false );
        integrationResultWithDissipationInIo = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

        for( std::map< double, Eigen::VectorXd >::const_iterator mapIterator = integrationResultWithDissipationInIo.begin( );
             mapIterator != integrationResultWithDissipationInIo.end( ); mapIterator++ )
        {
            integrationResultWithDissipationInIoKepler[ mapIterator->first ] =
                    orbital_element_conversions::convertCartesianToKeplerianElements< double >(
                        mapIterator->second, getBodyGravitationalParameter( "Jupiter" ) +
                        getBodyGravitationalParameter( "Io" ) );
            semiMajorAxesIoTide[ mapIterator->first ] = integrationResultWithDissipationInIoKepler[ mapIterator->first ]( 0 ) -
                    integrationResultWithDissipationInIoKepler.begin( )->second( 0 );
            eccentricitiesIoTide[ mapIterator->first ] = integrationResultWithDissipationInIoKepler[ mapIterator->first ]( 1 ) -
                    integrationResultWithDissipationInIoKepler.begin( )->second( 1 );

        }
    }



    Eigen::Vector6d intialKeplerElements = integrationResultWithDissipationInIoKepler.begin( )->second;
    double meanMotion = basic_astrodynamics::computeKeplerMeanMotion(
                intialKeplerElements( 0 ), getBodyGravitationalParameter( "Jupiter" ) + getBodyGravitationalParameter( "Io" ) );
    double orbitalPeriod = 2.0 * mathematical_constants::PI / meanMotion;
//    double jupiterForcingTime = 9.925 * 3600.0 * orbitalPeriod / ( 2.0 * std::fabs( orbitalPeriod - 9.925 * 3600.0 ) );
//    double jupiterQualityFactor = 1.0 / std::sin( 2.0 * mathematical_constants::PI * jupiterTimeLag / ( jupiterForcingTime ) );

    double ioQualityFactor = 1.0 / std::sin( 2.0 * mathematical_constants::PI * ( ioTimeLag / 2.0 ) / orbitalPeriod );

    std::vector< double > semiMajorAxisFitIoTide = linear_algebra::getLeastSquaresPolynomialFit(
                semiMajorAxesIoTide, boost::assign::list_of( 0 )( 1 ) );
    std::vector< double > eccentricityFitIoTide = linear_algebra::getLeastSquaresPolynomialFit(
                eccentricitiesIoTide, boost::assign::list_of( 0 )( 1 ) );

    double theoreticalSemiMajorAxisRateFromIoTide =  -21.0 * ioLoveNumber / ioQualityFactor * getBodyGravitationalParameter( "Jupiter" ) /
            getBodyGravitationalParameter( "Io" ) * std::pow( getAverageRadius( "Io" ) / intialKeplerElements( 0 ), 5.0 ) *
            intialKeplerElements( 0 ) * meanMotion * intialKeplerElements( 1 ) * intialKeplerElements( 1 );
    double theoreticaEccentricityRateFromIoTide =  - 21.0 / 2.0 * ioLoveNumber / ioQualityFactor * getBodyGravitationalParameter( "Jupiter" ) /
            getBodyGravitationalParameter( "Io" ) * std::pow( getAverageRadius( "Io" ) / intialKeplerElements( 0 ), 5.0 ) *
            intialKeplerElements( 1 ) * meanMotion;


    std::cout<<"THEORETICAL VALUE "<<theoreticalSemiMajorAxisRateFromIoTide<<" "<<theoreticaEccentricityRateFromIoTide<<" "<<std::endl<<
               ioQualityFactor<<" "<<getAverageRadius( "Jupiter" ) / intialKeplerElements( 0 )<<" "<<getBodyGravitationalParameter( "Io" ) /
               getBodyGravitationalParameter( "Jupiter" )<<std::endl;

    std::cout<<"COMPARISON "<<semiMajorAxisFitIoTide.at( 0 )<<" "<<semiMajorAxisFitIoTide.at( 1 )<<" "<<eccentricityFitIoTide.at( 0 )<<" "<<eccentricityFitIoTide.at( 1 )<<" "<<
               theoreticalSemiMajorAxisRateFromIoTide / semiMajorAxisFitIoTide.at( 1 )<<" "<<
               theoreticaEccentricityRateFromIoTide /  eccentricityFitIoTide.at( 1 )<<" "<<std::endl;



}

BOOST_AUTO_TEST_SUITE_END( )

}

}
