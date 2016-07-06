

#define BOOST_TEST_MAIN

#include <boost/assign/list_of.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/bind.hpp>

#include "Tudat/Astrodynamics/Relativity/relativisticLightTimeCorrection.h"
#include "Tudat/Astrodynamics/Ephemerides/constantEphemeris.h"
#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/firstOrderRelativisticLightTimeCorrection.h"

namespace tudat
{

namespace unit_tests
{

using namespace tudat::observation_models;
using namespace tudat::relativity;
using namespace tudat::ephemerides;

BOOST_AUTO_TEST_SUITE( test_shapiro_delay )

BOOST_AUTO_TEST_CASE( testShapiroDelay )
{
    basic_mathematics::Vector6d groundStationState;
    groundStationState<<0.0, 0.0, 6378.0, 0.0, 0.0, 0.0;
    basic_mathematics::Vector6d satelliteState;
    satelliteState << 0.0, 0.0, 26600.0, 0.0, 0.0, 0.0;
    basic_mathematics::Vector6d centralBodyPosition = basic_mathematics::Vector6d::Zero( );

    boost::shared_ptr< ConstantEphemeris > ephemeris = boost::make_shared< ConstantEphemeris >(
                boost::lambda::constant( centralBodyPosition ) );

    double earthGravitationalParameter = 398600.44189E9;

    double directCalculation = calculateFirstOrderLightTimeCorrectionFromCentralBody(
                earthGravitationalParameter, groundStationState.segment( 0, 3 ),
                satelliteState.segment( 0, 3 ), centralBodyPosition.segment( 0, 3 ) );

    std::vector< boost::function< basic_mathematics::Vector6d( const double ) > > perturbingBodyStateFunctions;
    std::vector< boost::function< double( ) > > perturbingBodyGravitationalParameterFunctions;

    perturbingBodyStateFunctions.push_back( boost::bind( &Ephemeris::getCartesianStateFromEphemeris, ephemeris, _1,
                                                         basic_astrodynamics::JULIAN_DAY_ON_J2000 ) );
    perturbingBodyGravitationalParameterFunctions.push_back( boost::lambda::constant( earthGravitationalParameter ) );

    FirstOrderLightTimeCorrectionCalculator correctionCalculator(
                perturbingBodyStateFunctions, perturbingBodyGravitationalParameterFunctions,
                boost::assign::list_of( "Earth" ), "Satellite", "Earth" );

    double classInterfaceCalculation = correctionCalculator.calculateLightTimeCorrection(
                groundStationState, satelliteState, 0.0, 0.0 );

    // Living reviews in relativity, GPS.
    double expectedResult = 6.3E-3;

    BOOST_CHECK_CLOSE_FRACTION( 0.5 * classInterfaceCalculation * physical_constants::SPEED_OF_LIGHT, expectedResult, 6.0E-2 );
    BOOST_CHECK_CLOSE_FRACTION( 0.5 * directCalculation * physical_constants::SPEED_OF_LIGHT, expectedResult, 6.0E-2 );
}

BOOST_AUTO_TEST_SUITE_END( )

}

}
