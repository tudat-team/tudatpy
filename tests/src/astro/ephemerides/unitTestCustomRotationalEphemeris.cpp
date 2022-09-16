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

#include <limits>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/basics/testMacros.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/environment_setup/createRotationModel.h"

namespace tudat
{
namespace unit_tests
{

using namespace ephemerides;
using namespace simulation_setup;
using namespace spice_interface;


BOOST_AUTO_TEST_SUITE( test_custom_rotational_ephemeris )

Eigen::Matrix3d getCustomRotationMatrix( const double time )
{
    return Eigen::Matrix3d( spice_interface::computeRotationQuaternionBetweenFrames(
                "IAU_Earth", "J2000", time ) );
}

// Test functions to calculate rotation matrix derivative from angular velocity vector and vice
// versa.
BOOST_AUTO_TEST_CASE( testCustomRotationalEphemeris )
{
    loadStandardSpiceKernels( );

    std::function< Eigen::Matrix3d( const double ) > customRotationFunction =
            &getCustomRotationMatrix;
    std::shared_ptr< RotationModelSettings > rotationSettings =
            std::make_shared< CustomRotationModelSettings >(
                "J2000", "IAU_Earth", customRotationFunction, 60.0 );

    std::shared_ptr< ephemerides::RotationalEphemeris > rotationModel = createRotationModel(
                rotationSettings, "Earth", SystemOfBodies( ) );

    double testTime = 1.0E7;
    Eigen::Matrix3d testRotation = rotationModel->getRotationMatrixToTargetFrame( testTime );
    Eigen::Matrix3d spiceRotation = Eigen::Matrix3d( computeRotationQuaternionBetweenFrames(
                                                         "J2000", "IAU_Earth", testTime ) );

    Eigen::Matrix3d testRotationDerivative = rotationModel->getDerivativeOfRotationToTargetFrame( testTime );
    Eigen::Matrix3d spiceRotationDerivative = computeRotationMatrixDerivativeBetweenFrames( "J2000", "IAU_Earth", testTime );
\
    Eigen::Vector3d testAngularVelocity = rotationModel->getRotationalVelocityVectorInBaseFrame( testTime );
    Eigen::Vector3d spiceAngularVelocity = getAngularVelocityVectorOfFrameInOriginalFrame( "J2000", "IAU_Earth", testTime );

    for( unsigned int i = 0; i < 3; i++ )
    {
        for( unsigned int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( ( testRotation( i, j ) - spiceRotation( i, j ) ), 2.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL( ( testRotationDerivative( i, j ) - spiceRotationDerivative( i, j ) ), 1.0E-9 );

        }
    }
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testAngularVelocity, spiceAngularVelocity, 1.0E-5 );

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
