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
 *      120821    K. Kumar          Unit test created.
 *      120827    K. Kumar          Cleaned up code and added reversed-order case.
 *      120913    K. Kumar          Rewrote unit test to use test state derivative models;
 *                                  implemented matrix-based and vector-based composite state
 *                                  tests.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <iostream>

#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/tuple/tuple.hpp>

#include <TudatCore/Basics/testMacros.h>

#include "Tudat/Astrodynamics/StateDerivativeModels/compositeStateDerivativeModel.h"
#include "Tudat/Astrodynamics/StateDerivativeModels/UnitTests/testStateDerivativeModels.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_composite_state_derivative_model )

//! Test whether composite state derivative model works correctly with matrices.
BOOST_AUTO_TEST_CASE( test_CompositeStateDerivativeModelWithMatrices )
{
    // Shortcuts.
    typedef state_derivative_models::CompositeStateDerivativeModel<
            double, Matrix2x3d, Eigen::MatrixXd > CompositeStateDerivativeModel;
    typedef boost::shared_ptr< CompositeStateDerivativeModel >
            CompositeStateDerivativeModelPointer;

    // This test creates and evalutes a composite state derivative model for system of coupled,
    // first-order, differential equations.

    // Set current time.
    time = 1.0;

    // Set current Vector2d state.
    vector2dState = Eigen::Vector2d( 1.3, -0.45 );

    // Set current Matrix2d state.
    matrix2dState = ( Eigen::MatrixXd( 2, 2 ) << 1.3, -0.45, 9.76, 4.56 ).finished( );

    // Set current composite state.
    Matrix2x3d currentCompositeState = ( Eigen::MatrixXd( 2, 3 )
                                         << vector2dState, matrix2dState ).finished( );

    // Create state derivative for matrix2d state.
    CoupledMatrix2dStateDerivativeModelPointer coupledMatrix2dStateDerivativeModel
            = boost::make_shared< CoupledMatrix2dStateDerivativeModel >( &getVector2dState );

    // Create state derivative model map ([Vector2d Matrix2d] composite state structure).
    CompositeStateDerivativeModel::StateDerivativeModelMap stateDerivativeModelMap;
    stateDerivativeModelMap[ boost::make_tuple( 0, 0, 2, 1 ) ] = &computeVector2dStateDerivative;
    stateDerivativeModelMap[ boost::make_tuple( 0, 1, 2, 2 ) ]
            = boost::bind( &CoupledMatrix2dStateDerivativeModel::computeStateDerivative,
                           coupledMatrix2dStateDerivativeModel, _1, _2 );

    // Create composite state derivative model.
    CompositeStateDerivativeModelPointer compositeStateDerivativeModel
            = boost::make_shared< CompositeStateDerivativeModel >(
                stateDerivativeModelMap, &updateMatrixData );

    // Compute composite state derivative.
    const Matrix2x3d computedCompositeStateDerivative
            = compositeStateDerivativeModel->computeStateDerivative( time, currentCompositeState );

    // Set expected composite state derivative.
    Matrix2x3d expectedCompositeStateDerivative;
    expectedCompositeStateDerivative.block( 0, 0, 2, 1 )
            = computeVector2dStateDerivative( time, vector2dState );
    expectedCompositeStateDerivative.block( 0, 1, 2, 2 )
            = coupledMatrix2dStateDerivativeModel->computeStateDerivative( time, matrix2dState );

    // Check that computed composite state derivative matches expected values.
    TUDAT_CHECK_MATRIX_BASE( computedCompositeStateDerivative,
                             expectedCompositeStateDerivative )
            BOOST_CHECK_EQUAL( computedCompositeStateDerivative.coeff( row, col ),
                               expectedCompositeStateDerivative.coeff( row, col ) );
}

//! Test whether composite state derivative model works correctly with vectors.
BOOST_AUTO_TEST_CASE( test_CompositeStateDerivativeModelWithVectors )
{
    // Shortcuts.
    typedef state_derivative_models::CompositeStateDerivativeModel<
            double, Vector5d, Eigen::VectorXd > CompositeStateDerivativeModel;
    typedef boost::shared_ptr< CompositeStateDerivativeModel >
            CompositeStateDerivativeModelPointer;

    // This test creates and evalutes a composite state derivative model for system of coupled,
    // first-order, differential equations.

    // Set current time.
    time = 2.2;

    // Set current Vector3d state.
    vector3dState = Eigen::Vector3d( 6.6, -1.13, 4.78 );

    // Set current Vector2d state.
    vector2dState = Eigen::Vector2d( 1.3, -0.45 );

    // Set current composite state.
    Vector5d currentCompositeState = ( Eigen::VectorXd( 5 )
                                       << vector3dState, vector2dState ).finished( );

    // Create state derivative model map ([Vector3d Vector2d] composite state structure).
    CompositeStateDerivativeModel::VectorStateDerivativeModelMap stateDerivativeModelMap;
    stateDerivativeModelMap[ std::make_pair( 0, 3 ) ] = &computeVector3dStateDerivative;
    stateDerivativeModelMap[ std::make_pair( 3, 2 ) ] = &computeVector2dStateDerivative;

    // Create composite state derivative model.
    CompositeStateDerivativeModelPointer compositeStateDerivativeModel
            = boost::make_shared< CompositeStateDerivativeModel >(
                stateDerivativeModelMap, &updateVectorData );

    // Compute composite state derivative.
    const Vector5d computedCompositeStateDerivative
            = compositeStateDerivativeModel->computeStateDerivative( time, currentCompositeState );

    // Set expected composite state derivative.
    Vector5d expectedCompositeStateDerivative;
    expectedCompositeStateDerivative.segment( 0, 3 )
            = computeVector3dStateDerivative( time, vector3dState );
    expectedCompositeStateDerivative.segment( 3, 2 )
            = computeVector2dStateDerivative( time, vector2dState );

    // Check that computed composite state derivative matches expected values.
    TUDAT_CHECK_MATRIX_BASE( computedCompositeStateDerivative,
                             expectedCompositeStateDerivative )
            BOOST_CHECK_EQUAL( computedCompositeStateDerivative.coeff( row, col ),
                               expectedCompositeStateDerivative.coeff( row, col ) );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
