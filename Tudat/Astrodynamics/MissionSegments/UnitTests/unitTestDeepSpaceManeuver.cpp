/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110527    K. Kumar          First creation of the code.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      120416    T. Secretin       Boostified unit test.
 *
 *    References
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <limits>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <TudatCore/Basics/testMacros.h>

#include "Tudat/Astrodynamics/MissionSegments/deepSpaceManeuver.h"
#include "Tudat/Astrodynamics/States/state.h"
#include "Tudat/Astrodynamics/States/cartesianElements.h"

namespace tudat
{
namespace unit_tests
{

//! Test deep space maneuver.
BOOST_AUTO_TEST_SUITE( test_deep_space_maneuver )

//! Test setTime( ) and getTime( ) functions.
BOOST_AUTO_TEST_CASE( testSetTimeAndGetTimeFunctions )
{
    // Declare and initialize time of DSM.
    const double timeOfDeepSpaceManeuver_ = 123.3;

    // Create DSM object and set time of DSM.
    using astrodynamics::mission_segments::DeepSpaceManeuver;
    DeepSpaceManeuver deepSpaceManeuver_( -0.0, timeOfDeepSpaceManeuver_,
                                          boost::shared_ptr< astrodynamics::states::State >( ) );

    BOOST_CHECK_CLOSE_FRACTION( timeOfDeepSpaceManeuver_, deepSpaceManeuver_.getTime( ),
                                std::numeric_limits< double >::epsilon( ) );
}

//! Test setState( ) and getState( ) functions with Cartesian elements state.
BOOST_AUTO_TEST_CASE ( testSetStateAndGetStateFunctionsWithCartesianElements )
{
    // Declare and initialize Cartesian elements state for DSM.
    using astrodynamics::states::CartesianElements;
    boost::shared_ptr< CartesianElements > deepSpaceManeuverState_
            = boost::make_shared< CartesianElements >( );
    deepSpaceManeuverState_->setCartesianElementX( 1000.3 );
    deepSpaceManeuverState_->setCartesianElementY( 12.3 );
    deepSpaceManeuverState_->setCartesianElementZ( 613.41 );
    deepSpaceManeuverState_->setCartesianElementXDot( 2300.5 );
    deepSpaceManeuverState_->setCartesianElementYDot( 100.34 );
    deepSpaceManeuverState_->setCartesianElementZ( 342.6 );

    // Create DSM object and set DSM state.
    using astrodynamics::mission_segments::DeepSpaceManeuver;
    DeepSpaceManeuver deepSpaceManeuver_( -0.0, -0.0, deepSpaceManeuverState_ );

    // Test if getState( ) function results in set state for DSM.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( deepSpaceManeuverState_->state,
                                           deepSpaceManeuver_.getState( )->state, 1.0e-15 );
}

//! Test setDeltaV( ) and getDeltaV( ) functions.
BOOST_AUTO_TEST_CASE ( testSetDeltaVandGetDeltaVFunctions )
{
    // Declare and initialize deltaV for DSM.
    const double deltaVOfDeepSpaceManeuver_ = 254.78;

    // Create DSM object and set deltaV of DSM.
    using astrodynamics::mission_segments::DeepSpaceManeuver;
    DeepSpaceManeuver deepSpaceManeuver_( deltaVOfDeepSpaceManeuver_, -0.0,
                                          boost::shared_ptr< astrodynamics::states::State >( ) );

    // Test if getDeltaV( ) function results in set deltaV for DSM.
    BOOST_CHECK_CLOSE_FRACTION( deltaVOfDeepSpaceManeuver_, deepSpaceManeuver_.getDeltaV( ),
                                std::numeric_limits< double >::epsilon( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
