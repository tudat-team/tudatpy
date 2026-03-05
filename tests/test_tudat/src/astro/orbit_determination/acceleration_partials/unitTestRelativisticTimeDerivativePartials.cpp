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

#include <functional>
#include <limits>
#include <string>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/astro/orbit_determination/acceleration_partials/numericalAccelerationPartial.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/gravitationalParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/ppnParameters.h"
#include "tudat/astro/propagators/relativisticTimeStateDerivative.h"
#include "tudat/astro/relativity/metric.h"
#include "tudat/astro/relativity/schwarzschildMetric.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createStateDerivativePartials.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::orbital_element_conversions;
using namespace tudat::unit_conversions;
using namespace tudat::estimatable_parameters;
using namespace tudat::relativity;
using namespace tudat::gravitation;
using namespace tudat::acceleration_partials;
using namespace tudat::basic_mathematics;

namespace
{

Eigen::MatrixXd getStatePartial(
        const std::pair< std::function< void( Eigen::Block< Eigen::MatrixXd > ) >, int >& partialData,
        const int numberOfRows )
{
    Eigen::MatrixXd partial = Eigen::MatrixXd::Zero( numberOfRows, partialData.second );
    if( partialData.second > 0 )
    {
        partialData.first( partial.block( 0, 0, numberOfRows, partialData.second ) );
    }
    return partial;
}

}  // namespace

BOOST_AUTO_TEST_SUITE( test_relativistic_state_derivative_partials )

BOOST_AUTO_TEST_CASE( testFirstOrderBarycentricToBodycentricPartials )
{
    spice_interface::loadStandardSpiceKernels( );

    const double evaluationTime = 1.0E6;

    auto sun = std::make_shared< Body >( );
    auto earth = std::make_shared< Body >( );
    auto moon = std::make_shared< Body >( );

    sun->setEphemeris(
        std::make_shared< ephemerides::ConstantEphemeris >( getBodyCartesianStateAtEpoch( "Sun", "SSB", "ECLIPJ2000", "NONE", evaluationTime ),
            "SSB", "ECLIPJ2000" ) );
    earth->setEphemeris(
    std::make_shared< ephemerides::ConstantEphemeris >( getBodyCartesianStateAtEpoch( "Earth", "SSB", "ECLIPJ2000", "NONE", evaluationTime ),
        "SSB", "ECLIPJ2000" ) );
    moon->setEphemeris(
    std::make_shared< ephemerides::ConstantEphemeris >( getBodyCartesianStateAtEpoch( "Moon", "SSB", "ECLIPJ2000", "NONE", evaluationTime ),
        "SSB", "ECLIPJ2000" ) );

    auto sunGravityFieldModel =
            std::make_shared< gravitation::GravityFieldModel >( spice_interface::getBodyGravitationalParameter( "Sun" ) );
    sun->setGravityFieldModel( sunGravityFieldModel );
    auto earthGravityFieldModel =
            std::make_shared< gravitation::GravityFieldModel >( spice_interface::getBodyGravitationalParameter( "Earth" ) );
    earth->setGravityFieldModel( earthGravityFieldModel );
    auto moonGravityFieldModel =
            std::make_shared< gravitation::GravityFieldModel >( spice_interface::getBodyGravitationalParameter( "Moon" ) );
    moon->setGravityFieldModel( moonGravityFieldModel );

    SystemOfBodies bodies;
    bodies.addBody( earth, "Earth" );
    bodies.addBody( moon, "Moon" );
    bodies.addBody( sun, "Sun" );

    std::vector< std::string > perturbingBodies{ "Sun", "Moon" };
    auto timeDerivativeModel = std::make_shared< FirstOrderBarycentricToBodyCentricTimeStateDerivative< double, double > >(
            bodies, "Earth", perturbingBodies );
    timeDerivativeModel->updateStateDerivativeModel( evaluationTime );

    auto timeDerivativePartial =
            createRelativisticTimeStateDerivativePartial< double, double, double >( timeDerivativeModel, bodies, nullptr );
    timeDerivativePartial->update( evaluationTime );

    Eigen::MatrixXd partialWrtEarthState = getStatePartial(
            timeDerivativePartial->getDerivativeFunctionWrtStateOfIntegratedBody(
                    std::make_pair( "Earth", "" ), translational_state ),
            1 );
    Eigen::MatrixXd partialWrtSunState = getStatePartial(
            timeDerivativePartial->getDerivativeFunctionWrtStateOfIntegratedBody(
                    std::make_pair( "Sun", "" ), translational_state ),
            1 );
    Eigen::MatrixXd partialWrtMoonState = getStatePartial(
            timeDerivativePartial->getDerivativeFunctionWrtStateOfIntegratedBody(
                    std::make_pair( "Moon", "" ), translational_state ),
            1 );

    auto sunGravitationalParameter = std::make_shared< GravitationalParameter >( sunGravityFieldModel, "Sun" );
    auto moonGravitationalParameter = std::make_shared< GravitationalParameter >( moonGravityFieldModel, "Moon" );

    Eigen::MatrixXd partialWrtSunGravitationalParameter = timeDerivativePartial->wrtParameter( sunGravitationalParameter );
    Eigen::MatrixXd partialWrtMoonGravitationalParameter = timeDerivativePartial->wrtParameter( moonGravitationalParameter );

    const Eigen::Vector6d statePerturbation = ( Eigen::Vector6d( ) << 100000.0, 100000.0, 100000.0, 1.0, 1.0, 1.0 ).finished( );

    const std::function< void( const Eigen::Vector6d& ) > earthStateSetFunction =
            [earth]( const Eigen::Vector6d& state ){ earth->setState( state ); };
    const std::function< void( const Eigen::Vector6d& ) > sunStateSetFunction =
            [sun]( const Eigen::Vector6d& state ){ sun->setState( state ); };
    const std::function< void( const Eigen::Vector6d& ) > moonStateSetFunction =
            [moon]( const Eigen::Vector6d& state ){ moon->setState( state ); };

    Eigen::Matrix< double, 1, 6 > testPartialWrtEarthState = calculateRelativisticTimeDerivativeWrtStatePartials(
            earthStateSetFunction, timeDerivativeModel, earth->getState( ), statePerturbation );
    Eigen::Matrix< double, 1, 6 > testPartialWrtSunState = calculateRelativisticTimeDerivativeWrtStatePartials(
            sunStateSetFunction, timeDerivativeModel, sun->getState( ), statePerturbation );
    Eigen::Matrix< double, 1, 6 > testPartialWrtMoonState = calculateRelativisticTimeDerivativeWrtStatePartials(
            moonStateSetFunction, timeDerivativeModel, moon->getState( ), statePerturbation );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthState, partialWrtEarthState, 2.0e-4 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunState, partialWrtSunState, 2.0e-4 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMoonState, partialWrtMoonState, 2.0e-4 );

    Eigen::Matrix< double, 1, Eigen::Dynamic > testPartialWrtSunGravitationalParameter =
            calculateRelativisticTimeDerivativeWrtParameterPartials(
                    sunGravitationalParameter, timeDerivativeModel, 1.0E-6 * sunGravitationalParameter->getParameterValue( ) );
    Eigen::Matrix< double, 1, Eigen::Dynamic > testPartialWrtMoonGravitationalParameter =
            calculateRelativisticTimeDerivativeWrtParameterPartials(
                    moonGravitationalParameter, timeDerivativeModel, 1.0E-6 * moonGravitationalParameter->getParameterValue( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            testPartialWrtSunGravitationalParameter, partialWrtSunGravitationalParameter, 1.0e-8 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            testPartialWrtMoonGravitationalParameter, partialWrtMoonGravitationalParameter, 1.0e-5 );
}

BOOST_AUTO_TEST_CASE( testDirectProperTimeDerivativePartials )
{
    spice_interface::loadStandardSpiceKernels( );

    auto sun = std::make_shared< Body >( );
    sun->setState( getBodyCartesianStateAtEpoch( "Sun", "SSB", "ECLIPJ2000", "NONE", 1.0E6 ) );

    auto sunGravityFieldModel =
            std::make_shared< gravitation::GravityFieldModel >( spice_interface::getBodyGravitationalParameter( "Sun" ) );
    sun->setGravityFieldModel( sunGravityFieldModel );

    auto testParticle = std::make_shared< Body >( );
    Eigen::Vector6d initialKeplerElements;
    initialKeplerElements[ semiMajorAxisIndex ] = 50.0E9;
    initialKeplerElements[ eccentricityIndex ] = 0.5;
    initialKeplerElements[ inclinationIndex ] = 30.0 * mathematical_constants::PI / 180.0;
    initialKeplerElements[ argumentOfPeriapsisIndex ] = 1.0;
    initialKeplerElements[ longitudeOfAscendingNodeIndex ] = 1.0;
    initialKeplerElements[ trueAnomalyIndex ] = 0.2;
    testParticle->setState( convertKeplerianToCartesianElements(
            initialKeplerElements, sunGravityFieldModel->getGravitationalParameter( ) ) );

    SystemOfBodies bodies;
    bodies.addBody( sun, "Sun" );
    bodies.addBody( testParticle, "TestParticle" );

    auto ppnParameterSet = std::make_shared< PPNParameterSet >( 1.0, 1.0 );
    auto baseMetric = std::make_shared< HarmonicSchwarzschildMetric >(
            std::bind( &GravityFieldModel::getGravitationalParameter, sunGravityFieldModel ),
            ppnParameterSet,
            "Sun" );

    auto timeDerivativeModel = std::make_shared< DirectProperTimeRateStateDerivative< double, double > >(
            baseMetric, std::make_pair( "TestParticle", "" ), bodies );
    timeDerivativeModel->updateStateDerivativeModel( 0.0 );

    auto timeDerivativePartial =
            createRelativisticTimeStateDerivativePartial< double, double, double >( timeDerivativeModel, bodies, nullptr );
    timeDerivativePartial->update( 0.0 );

    Eigen::MatrixXd partialWrtTestParticleState = getStatePartial(
            timeDerivativePartial->getDerivativeFunctionWrtStateOfIntegratedBody(
                    std::make_pair( "TestParticle", "" ), translational_state ),
            1 );

    auto sunGravitationalParameter = std::make_shared< GravitationalParameter >( sunGravityFieldModel, "Sun" );
    auto ppnParameterGamma = std::make_shared< PPNParameterGamma >( ppnParameterSet );
    auto ppnParameterBeta = std::make_shared< PPNParameterBeta >( ppnParameterSet );

    Eigen::MatrixXd partialWrtSunGravitationalParameter = timeDerivativePartial->wrtParameter( sunGravitationalParameter );
    Eigen::MatrixXd partialWrtPpnGamma = timeDerivativePartial->wrtParameter( ppnParameterGamma );
    Eigen::MatrixXd partialWrtPpnBeta = timeDerivativePartial->wrtParameter( ppnParameterBeta );

    const Eigen::Vector6d statePerturbation = ( Eigen::Vector6d( ) << 100000.0, 100000.0, 100000.0, 1.0, 1.0, 1.0 ).finished( );

    const std::function< void( const Eigen::Vector6d& ) > testParticleStateSetFunction =
            [testParticle]( const Eigen::Vector6d& state ){ testParticle->setState( state ); };

    Eigen::Matrix< double, 1, 6 > testPartialWrtTestParticleState = calculateRelativisticTimeDerivativeWrtStatePartials(
            testParticleStateSetFunction, timeDerivativeModel, testParticle->getState( ), statePerturbation );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( -1.0 * testPartialWrtTestParticleState ), partialWrtTestParticleState, 1.0e-9 );

    Eigen::Matrix< double, 1, Eigen::Dynamic > testPartialWrtSunGravitationalParameter =
            calculateRelativisticTimeDerivativeWrtParameterPartials(
                    sunGravitationalParameter, timeDerivativeModel, 1.0E-6 * sunGravitationalParameter->getParameterValue( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            ( -1.0 * testPartialWrtSunGravitationalParameter ), partialWrtSunGravitationalParameter, 1.0e-9 );

    Eigen::Matrix< double, 1, Eigen::Dynamic > testPartialWrtPpnGamma =
            calculateRelativisticTimeDerivativeWrtParameterPartials( ppnParameterGamma, timeDerivativeModel, 0.1 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( -1.0 * testPartialWrtPpnGamma ), partialWrtPpnGamma, 1.0e-6 );

    Eigen::Matrix< double, 1, Eigen::Dynamic > testPartialWrtPpnBeta =
            calculateRelativisticTimeDerivativeWrtParameterPartials( ppnParameterBeta, timeDerivativeModel, 0.1 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( -1.0 * testPartialWrtPpnBeta ), partialWrtPpnBeta, 1.0e-6 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
