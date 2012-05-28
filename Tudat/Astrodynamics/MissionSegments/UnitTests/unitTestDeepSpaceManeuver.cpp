/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
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
