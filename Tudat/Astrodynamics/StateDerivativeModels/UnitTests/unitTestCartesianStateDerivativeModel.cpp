/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      120820    K. Kumar          Unit test created.
 *      120911    K. Kumar          Added unit test to demonstrate frame transformation
 *                                  functionality.
 *
 *    References
 *
 *    Notes:
 *      Test tolerance was set at 5.0e-15 instead of epsilon due to rounding errors in Eigen types
 *      with entries over a number of orders of magnitude, presumably causing the observed larger
 *      than epsilon relative differences between expected and computed values.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/assign/list_of.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/UnitTests/testAccelerationModels.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/UnitTests/testBody.h"
#include "Tudat/Astrodynamics/StateDerivativeModels/cartesianStateDerivativeModel.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

namespace tudat
{
namespace unit_tests
{

using boost::assign::list_of;
using Eigen::Vector6d;

//! Rotate vector over arbitrary angles.
/*!
 * This function computes a composite rotation of the input vector over arbitrary angles about the
 * unit x-, y-, and z-axes. This is used in conjunction with rotateOverOtherArbitraryAngles() to
 * test the frame transformation features provided by the CartesianStateDerivativeModel.
 * \param inputVector Input vector before rotation.
 * \return Rotated vector.
 * \sa CartesianStateDerivativeModel, rotateOverOtherArbitraryAngles().
 */
Eigen::Vector3d rotateOverArbitraryAngles( const Eigen::Vector3d& inputVector )
{
    // Declare rotation matrix.
    Eigen::Matrix3d rotationMatrix;
    rotationMatrix = Eigen::AngleAxisd( -1.15, Eigen::Vector3d::UnitX( ) )
            * Eigen::AngleAxisd( 0.23, Eigen::Vector3d::UnitY( ) )
            * Eigen::AngleAxisd( 2.56, Eigen::Vector3d::UnitZ( ) );

    // Compute rotated matrix and return.
    return  rotationMatrix * inputVector;
}

//! Rotate vector over other arbitrary angles.
/*!
 * This function computes another composite rotation of the input vector over different arbitrary
 * angles (compared to the rotateOverArbitraryAngles() function) about the unit x-, y-, and z-axes.
 * This is used in conjunction with rotateOverArbitraryAngles() to test the frame transformation
 * features provided by the CartesianStateDerivativeModel.
 * \param inputVector Input vector before rotation.
 * \return Rotated vector.
 * \sa CartesianStateDerivativeModel, rotateOverArbitraryAngles().
 */
Eigen::Vector3d rotateOverOtherArbitraryAngles( const Eigen::Vector3d& inputVector )
{
    // Declare rotation matrix.
    Eigen::Matrix3d rotationMatrix;
    rotationMatrix = Eigen::AngleAxisd( 0.24, Eigen::Vector3d::UnitX( ) )
            * Eigen::AngleAxisd( -1.55, Eigen::Vector3d::UnitY( ) )
            * Eigen::AngleAxisd( 2.13, Eigen::Vector3d::UnitZ( ) );

    // Compute rotated matrix and return.
    return  rotationMatrix * inputVector;
}

BOOST_AUTO_TEST_SUITE( test_cartesian_state_derivative_model )

//! Test whether 6D Cartesian state derivative model works correctly without frame transformations.
BOOST_AUTO_TEST_CASE( test_CartesianStateDerivativeModel6DWithoutFrameTransformations )
{
    using basic_astrodynamics::AccelerationModel3dPointer;
    using state_derivative_models::CartesianStateDerivativeModel6d;
    using state_derivative_models::CartesianStateDerivativeModel6dPointer;

    // Shortcuts.
    typedef TestBody< 3, double > TestBody3d;
    typedef boost::shared_ptr< TestBody3d > TestBody3dPointer;
    typedef DerivedAccelerationModel< > DerivedAccelerationModel3d;
    typedef AnotherDerivedAccelerationModel< > AnotherDerivedAccelerationModel3d;

    // Set current state.
    const Vector6d currentState = ( Eigen::Vector6d( )
                                    << Eigen::Vector3d( -1.1, 2.2, -3.3 ),
                                    Eigen::Vector3d( 0.23, 1.67, -0.11 ) ).finished( );

    // Set current time.
    const double currentTime = 5.6;

    // Set current position.
    const Eigen::Vector3d currentPosition = currentState.segment( 0, 3 );

    // Set current velocity.
    const Eigen::Vector3d currentVelocity = currentState.segment( 3, 3 );

    // Create body with zombie time and state.
    TestBody3dPointer body = boost::make_shared< TestBody3d >( Eigen::VectorXd::Zero( 6 ), 0.0 );

    // Create acceleration models.
    AccelerationModel3dPointer firstAccelerationModel3d
            = boost::make_shared< DerivedAccelerationModel3d >(
                boost::bind( &TestBody3d::getCurrentPosition, body ),
                boost::bind( &TestBody3d::getCurrentTime, body ) );

    AccelerationModel3dPointer secondAccelerationModel3d
            = boost::make_shared< AnotherDerivedAccelerationModel3d >(
                boost::bind( &TestBody3d::getCurrentPosition, body ),
                boost::bind( &TestBody3d::getCurrentVelocity, body ),
                boost::bind( &TestBody3d::getCurrentTime, body ) );

    // Create list of acceleration models to provide to state derivative model.
    CartesianStateDerivativeModel6d::AccelerationModelPointerVector listOfAccelerations
            = list_of( firstAccelerationModel3d )( secondAccelerationModel3d );

    // Declare Cartesian state derivative model.
    CartesianStateDerivativeModel6dPointer stateDerivativeModel
            = boost::make_shared< CartesianStateDerivativeModel6d >(
                listOfAccelerations,
                boost::bind( &TestBody3d::setCurrentTimeAndState, body, _1, _2 ) );

    // Set expected accelerations.
    const Eigen::Vector3d expectedAccelerationFirstModel
            = currentPosition / ( currentTime * currentTime );
    const Eigen::Vector3d expectedAccelerationSecondModel
            = 0.5 * currentPosition / ( 3.2 * ( currentTime + 3.4 ) * currentTime )
            + currentVelocity / currentTime;

    // Set expected (cumulative) Cartesian state derivative.
    const Vector6d expectedCartesianStateDerivative
            = ( Eigen::Vector6d( ) << currentVelocity,
                expectedAccelerationFirstModel + expectedAccelerationSecondModel ).finished( );

    // Compute Cartesian state derivative.
    const Vector6d computedCartesianStateDerivative
            = stateDerivativeModel->computeStateDerivative( currentTime, currentState );

    // Check that computed Cartesian state derivative matches expected values.
    TUDAT_CHECK_MATRIX_BASE( computedCartesianStateDerivative, expectedCartesianStateDerivative )
            BOOST_CHECK_SMALL( computedCartesianStateDerivative.coeff( row, col ) -
                               expectedCartesianStateDerivative.coeff( row, col ),
                               5.0e-15 );
}

//! Test whether 6D Cartesian state derivative model works correctly with frame transformations.
BOOST_AUTO_TEST_CASE( test_CartesianStateDerivativeModel6DWithFrameTransformations )
{
    using basic_astrodynamics::AccelerationModel3dPointer;
    using state_derivative_models::CartesianStateDerivativeModel6d;
    using state_derivative_models::CartesianStateDerivativeModel6dPointer;

    // Shortcuts.
    typedef TestBody< 3, double > TestBody3d;
    typedef boost::shared_ptr< TestBody3d > TestBody3dPointer;
    typedef DerivedAccelerationModel< > DerivedAccelerationModel3d;
    typedef AnotherDerivedAccelerationModel< > AnotherDerivedAccelerationModel3d;

    // Set current state.
    const Vector6d currentState = ( Eigen::Vector6d( )
                                    << Eigen::Vector3d( -1.1, 2.2, -3.3 ),
                                    Eigen::Vector3d( 0.23, 1.67, -0.11 ) ).finished( );

    // Set current time.
    const double currentTime = 5.6;

    // Set current position.
    const Eigen::Vector3d currentPosition = currentState.segment( 0, 3 );

    // Set current velocity.
    const Eigen::Vector3d currentVelocity = currentState.segment( 3, 3 );

    // Create body with zombie time and state.
    TestBody3dPointer body = boost::make_shared< TestBody3d >( Eigen::VectorXd::Zero( 6 ), 0.0 );

    // Create acceleration models.
    AccelerationModel3dPointer firstAccelerationModel3d
            = boost::make_shared< DerivedAccelerationModel3d >(
                boost::bind( &TestBody3d::getCurrentPosition, body ),
                boost::bind( &TestBody3d::getCurrentTime, body ) );

    AccelerationModel3dPointer secondAccelerationModel3d
            = boost::make_shared< AnotherDerivedAccelerationModel3d >(
                boost::bind( &TestBody3d::getCurrentPosition, body ),
                boost::bind( &TestBody3d::getCurrentVelocity, body ),
                boost::bind( &TestBody3d::getCurrentTime, body ) );

    // Create list of reference frame transformations for first acceleration model.
    // NB: The order of this list is VERY important! The order of transformations executed is from
    // the beginning of the vector to the end sequentially.
    CartesianStateDerivativeModel6d::ListOfReferenceFrameTransformations listOfFrameTransformations
            = list_of( &rotateOverArbitraryAngles )( &rotateOverOtherArbitraryAngles );

    // Create list to pass to constructor of acceleration model/frame transformation list pairs.
    // In this case, there are two frame transformations executed for the first acceleration model,
    // and none for the second.
    CartesianStateDerivativeModel6d::ListOfAccelerationFrameTransformationPairs
            listOfAccelerationFrameTransformations
            = list_of( std::make_pair( firstAccelerationModel3d, listOfFrameTransformations ) )
            ( std::make_pair( secondAccelerationModel3d,
                              CartesianStateDerivativeModel6d::
                              ListOfReferenceFrameTransformations( ) ) );

    // Declare Cartesian state derivative model.
    CartesianStateDerivativeModel6dPointer stateDerivativeModel
            = boost::make_shared< CartesianStateDerivativeModel6d >(
                listOfAccelerationFrameTransformations,
                boost::bind( &TestBody3d::setCurrentTimeAndState, body, _1, _2 ) );

    // Set expected accelerations.
    const Eigen::Vector3d expectedAccelerationFirstModel
            = rotateOverOtherArbitraryAngles(
                rotateOverArbitraryAngles( currentPosition / ( currentTime * currentTime ) ) );
    const Eigen::Vector3d expectedAccelerationSecondModel
            = 0.5 * currentPosition / ( 3.2 * ( currentTime + 3.4 ) * currentTime )
            + currentVelocity / currentTime;

    // Set expected (cumulative) Cartesian state derivative.
    const Vector6d expectedCartesianStateDerivative
            = ( Eigen::Vector6d( ) << currentVelocity,
                expectedAccelerationFirstModel + expectedAccelerationSecondModel ).finished( );

    // Compute Cartesian state derivative.
    const Vector6d computedCartesianStateDerivative
            = stateDerivativeModel->computeStateDerivative( currentTime, currentState );

    // Check that computed Cartesian state derivative matches expected values.
    TUDAT_CHECK_MATRIX_BASE( computedCartesianStateDerivative, expectedCartesianStateDerivative )
            BOOST_CHECK_CLOSE_FRACTION( computedCartesianStateDerivative.coeff( row, col ),
                                        expectedCartesianStateDerivative.coeff( row, col ),
                                        5.0e-15 );
}

//! Test whether 4D Cartesian state derivative model works correctly.
BOOST_AUTO_TEST_CASE( test_CartesianStateDerivativeModel4D )
{
    using basic_astrodynamics::AccelerationModel2dPointer;
    using state_derivative_models::CartesianStateDerivativeModel4d;
    using state_derivative_models::CartesianStateDerivativeModel4dPointer;

    // Shortcuts.
    typedef TestBody< 2, double > TestBody2d;
    typedef boost::shared_ptr< TestBody2d > TestBody2dPointer;
    typedef DerivedAccelerationModel< Eigen::Vector2d, Eigen::Vector2d, double >
            DerivedAccelerationModel2d;
    typedef AnotherDerivedAccelerationModel< Eigen::Vector2d, Eigen::Vector2d,
            Eigen::Vector2d, double > AnotherDerivedAccelerationModel2d;

    // Set current state.
    const Eigen::Vector4d currentState( 3.45, 0.98, -0.12, -1.1e-4 );

    // Set current time.
    const double currentTime = -1.3;

    // Set current position.
    const Eigen::Vector2d currentPosition = currentState.segment( 0, 2 );

    // Set current velocity.
    const Eigen::Vector2d currentVelocity = currentState.segment( 2, 2 );

    // Create body with zombie time and state.
    TestBody2dPointer body = boost::make_shared< TestBody2d >( Eigen::VectorXd::Zero( 4 ), 0.0 );

    // Create acceleration model.
    AccelerationModel2dPointer firstAccelerationModel2d
            = boost::make_shared< DerivedAccelerationModel2d >(
                boost::bind( &TestBody2d::getCurrentPosition, body ),
                boost::bind( &TestBody2d::getCurrentTime, body ) );

    AccelerationModel2dPointer secondAccelerationModel2d
            = boost::make_shared< AnotherDerivedAccelerationModel2d >(
                boost::bind( &TestBody2d::getCurrentPosition, body ),
                boost::bind( &TestBody2d::getCurrentVelocity, body ),
                boost::bind( &TestBody2d::getCurrentTime, body ) );

    // Create list of acceleration models to provide to state derivative model.
    CartesianStateDerivativeModel4d::AccelerationModelPointerVector listOfAccelerations
            = list_of( firstAccelerationModel2d )( secondAccelerationModel2d );

    // Declare Cartesian state derivative model.
    CartesianStateDerivativeModel4dPointer stateDerivativeModel
            = boost::make_shared< CartesianStateDerivativeModel4d >(
                listOfAccelerations,
                boost::bind( &TestBody2d::setCurrentTimeAndState, body, _1, _2 ) );

    // Set expected accelerations.
    const Eigen::Vector2d expectedAccelerationFirstModel
            = currentPosition / ( currentTime * currentTime );
    const Eigen::Vector2d expectedAccelerationSecondModel
            = 0.5 * currentPosition / ( 3.2 * ( currentTime + 3.4 ) * currentTime )
            + currentVelocity / currentTime;

    // Set expected (cumulative) Cartesian state derivative.
    const Eigen::Vector4d expectedCartesianStateDerivative
            = ( Eigen::VectorXd( 4 ) << currentVelocity,
                expectedAccelerationFirstModel + expectedAccelerationSecondModel ).finished( );

    // Compute Cartesian state derivative.
    const Eigen::Vector4d computedCartesianStateDerivative
            = stateDerivativeModel->computeStateDerivative( currentTime, currentState );

    // Check that computed Cartesian state derivative matches expected values.
    TUDAT_CHECK_MATRIX_BASE( computedCartesianStateDerivative, expectedCartesianStateDerivative )
            BOOST_CHECK_CLOSE_FRACTION( computedCartesianStateDerivative.coeff( row, col ),
                                        expectedCartesianStateDerivative.coeff( row, col ),
                                        5.0e-15 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
