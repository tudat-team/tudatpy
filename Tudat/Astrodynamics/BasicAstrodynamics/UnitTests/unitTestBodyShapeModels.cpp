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
 *      150409    D. Dirkx          Migrated from personal code.
 *
 *    References
 *      Montebruck O, Gill E. Satellite Orbits, Springer, 2000.
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>

#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Basics/testMacros.h>

#include "Tudat/Astrodynamics/BasicAstrodynamics/oblateSpheroidBodyShapeModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalBodyShapeModel.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_geodetic_coordinate_conversions )

//! Test shape models, for test data, see testGeodeticCoordinateConversions.
BOOST_AUTO_TEST_CASE( testShapeModels )
{
    using namespace tudat::physical_constants;
    using namespace tudat::coordinate_conversions;
    using namespace tudat::unit_conversions;
    using namespace tudat::basic_astrodynamics;

    // Expected Cartesian state, Montenbruck & Gill (2000) Exercise 5.3.
    const Eigen::Vector3d testCartesianPosition( 1917032.190, 6029782.349, -801376.113 );

    // Expected Cartesian state, Montenbruck & Gill (2000) Exercise 5.3.
    const Eigen::Vector3d testGeodeticPosition( -63.667,
                                                convertDegreesToRadians( -7.26654999 ),
                                                convertDegreesToRadians( 72.36312094 ) );

    // Central body characteristics (WGS84 Earth ellipsoid).
    const double flattening = 1.0 / 298.257223563;
    const double equatorialRadius = 6378137.0;

    // Test spherealtitude
    {
        SphericalBodyShapeModel shapeModel = SphericalBodyShapeModel(
                    equatorialRadius );

        // Calculate altitude from shape object.
        const double altitudeFromObject = shapeModel.getAltitude(
                    testCartesianPosition );

        // Calculate object directly.
        const double directAltitude = testCartesianPosition.norm( ) - equatorialRadius;

        // Compare values.
        BOOST_CHECK_EQUAL( altitudeFromObject, directAltitude );
    }

    // Test oblate spheroid
    {
        OblateSpheroidBodyShapeModel shapeModel = OblateSpheroidBodyShapeModel(
                    equatorialRadius, flattening );

        // Calculate altitude from shape object.
        const double altitudeFromObject = shapeModel.getAltitude(
                    testCartesianPosition );

        // Calculate object from free function.
        const double directAltitude = calculateAltitudeOverOblateSpheroid(
                    testCartesianPosition, equatorialRadius, flattening, 1.0E-4 );

        // Compare values.
        BOOST_CHECK_SMALL( altitudeFromObject - testGeodeticPosition.x( ), 1.0E-4 );
        BOOST_CHECK_EQUAL( altitudeFromObject, directAltitude );

        // Test calculation of full geodetic position.
        Eigen::Vector3d calculatedGeodeticPosition = shapeModel.getGeodeticPositionWrtShape(
                    testCartesianPosition );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    calculatedGeodeticPosition, testGeodeticPosition, 1.0E-6 );
    }

    // Test free function altitude calculations
    {
        boost::shared_ptr< OblateSpheroidBodyShapeModel > shapeModel =
                boost::make_shared< OblateSpheroidBodyShapeModel >(
                    equatorialRadius, flattening );

        Eigen::Vector3d bodyPosition =
                ( Eigen::Vector3d( )<<
                  ASTRONOMICAL_UNIT / sqrt( 2.0 ),
                  ASTRONOMICAL_UNIT / sqrt( 2.0 ),
                  ASTRONOMICAL_UNIT  * 0.01 ).finished( );
        Eigen::Vector3d inertialTestCartesianPosition = testCartesianPosition + bodyPosition;

        double calculatedAltitute = getAltitudeFromNonBodyFixedPosition(
                    shapeModel, testCartesianPosition, Eigen::Vector3d::Zero( ),
                    Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) );
        BOOST_CHECK_SMALL( calculatedAltitute - testGeodeticPosition.x( ), 1.0E-4 );

        calculatedAltitute = getAltitudeFromNonBodyFixedPosition(
                    shapeModel, inertialTestCartesianPosition, bodyPosition,
                    Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) );
        BOOST_CHECK_SMALL( calculatedAltitute - testGeodeticPosition.x( ), 1.0E-4 );

        Eigen::Quaterniond dummyTestRotation =
                Eigen::AngleAxisd( 0.4343, Eigen::Vector3d::UnitZ( ) ) *
                Eigen::AngleAxisd( 2.4354, Eigen::Vector3d::UnitX( ) ) *
                Eigen::AngleAxisd( 1.2434, Eigen::Vector3d::UnitY( ) );

        calculatedAltitute = getAltitudeFromNonBodyFixedPosition(
                    shapeModel, dummyTestRotation * inertialTestCartesianPosition,
                    dummyTestRotation * bodyPosition,
                    dummyTestRotation.inverse( ) );
        BOOST_CHECK_SMALL( calculatedAltitute - testGeodeticPosition.x( ), 1.0E-4 );

        calculatedAltitute = getAltitudeFromNonBodyFixedPositionFunctions(
                    shapeModel, dummyTestRotation * inertialTestCartesianPosition,
                    boost::lambda::constant(
                        Eigen::Vector3d( dummyTestRotation * bodyPosition ) ),
                    boost::lambda::constant( dummyTestRotation.inverse( ) ) );
        BOOST_CHECK_SMALL( calculatedAltitute - testGeodeticPosition.x( ), 1.0E-4 );

    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

