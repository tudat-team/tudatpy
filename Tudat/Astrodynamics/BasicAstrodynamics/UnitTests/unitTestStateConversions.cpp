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
 *      130301    D. Dirkx          Migrated from personal code.
 *      130308    E.D. Brandon      Minor changes.
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

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/oblateSpheroidBodyShapeModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateRepresentationConversions.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_position_representation_conversions )

BOOST_AUTO_TEST_CASE( testGeneralGeodeticCoordinateConversions )
{
    using namespace coordinate_conversions;
    using namespace unit_conversions;

    // Expected Cartesian state, Montenbruck & Gill (2000) Exercise 5.3.
    const Eigen::Vector3d testCartesianPosition( 1917032.190, 6029782.349, -801376.113 );

    // Expected Cartesian state, Montenbruck & Gill (2000) Exercise 5.3.
    const Eigen::Vector3d testGeodeticPosition( -63.667,
                                                convertDegreesToRadians( -7.26654999 ),
                                                convertDegreesToRadians( 72.36312094 ) );

    Eigen::Vector3d testSphericalPosition = coordinate_conversions::convertCartesianToSpherical(
                testCartesianPosition );
    testSphericalPosition( 1 ) = mathematical_constants::PI / 2.0 - testSphericalPosition( 1 );


    // Central body characteristics (WGS84 Earth ellipsoid).
    const double flattening = 1.0 / 298.257223563;
    const double equatorialRadius = 6378137.0;

    boost::shared_ptr< basic_astrodynamics::OblateSpheroidBodyShapeModel > oblateSpheroidModel =
            boost::make_shared< basic_astrodynamics::OblateSpheroidBodyShapeModel >(
                equatorialRadius, flattening );

    Eigen::Vector3d computedCartesianPosition, computedSphericalPosition, computedGeodeticPosition;

    computedCartesianPosition = coordinate_conversions::convertPositionElements(
                testGeodeticPosition, coordinate_conversions::geodetic_position, coordinate_conversions::cartesian_position,
                oblateSpheroidModel );

    computedSphericalPosition = coordinate_conversions::convertPositionElements(
                testGeodeticPosition, coordinate_conversions::geodetic_position, coordinate_conversions::spherical_position,
                oblateSpheroidModel );


    for( unsigned int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( testCartesianPosition( i ) - computedCartesianPosition( i ), 1.0E-3 );
    }

    BOOST_CHECK_SMALL( testSphericalPosition( 0 ) - computedSphericalPosition( 0 ), 1.0E-3 );
    BOOST_CHECK_SMALL( testSphericalPosition( 1 ) - computedSphericalPosition( 1 ), 1.0E-3 / equatorialRadius );
    BOOST_CHECK_SMALL( testSphericalPosition( 2 ) - computedSphericalPosition( 2 ), 1.0E-3 / equatorialRadius );



    computedSphericalPosition = coordinate_conversions::convertPositionElements(
                testCartesianPosition, coordinate_conversions::cartesian_position, coordinate_conversions::spherical_position,
                oblateSpheroidModel );

    computedGeodeticPosition = coordinate_conversions::convertPositionElements(
                testCartesianPosition, coordinate_conversions::cartesian_position, coordinate_conversions::geodetic_position,
                oblateSpheroidModel );

    BOOST_CHECK_SMALL( testGeodeticPosition( 0 ) - computedGeodeticPosition( 0 ), 1.0E-3 );
    BOOST_CHECK_SMALL( testGeodeticPosition( 1 ) - computedGeodeticPosition( 1 ), 1.0E-3 / equatorialRadius );
    BOOST_CHECK_SMALL( testGeodeticPosition( 2 ) - computedGeodeticPosition( 2 ), 1.0E-3 / equatorialRadius );

    BOOST_CHECK_SMALL( testSphericalPosition( 0 ) - computedSphericalPosition( 0 ), 1.0E-3 );
    BOOST_CHECK_SMALL( testSphericalPosition( 1 ) - computedSphericalPosition( 1 ), 1.0E-3 / equatorialRadius );
    BOOST_CHECK_SMALL( testSphericalPosition( 2 ) - computedSphericalPosition( 2 ), 1.0E-3 / equatorialRadius );

    computedCartesianPosition = coordinate_conversions::convertPositionElements(
                testSphericalPosition, coordinate_conversions::spherical_position, coordinate_conversions::cartesian_position,
                oblateSpheroidModel );

    computedGeodeticPosition = coordinate_conversions::convertPositionElements(
                testSphericalPosition, coordinate_conversions::spherical_position, coordinate_conversions::geodetic_position,
                oblateSpheroidModel );

    for( unsigned int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( testCartesianPosition( i ) - computedCartesianPosition( i ), 1.0E-3 );
    }

    BOOST_CHECK_SMALL( testGeodeticPosition( 0 ) - computedGeodeticPosition( 0 ), 1.0E-3 );
    BOOST_CHECK_SMALL( testGeodeticPosition( 1 ) - computedGeodeticPosition( 1 ), 1.0E-3 / equatorialRadius );
    BOOST_CHECK_SMALL( testGeodeticPosition( 2 ) - computedGeodeticPosition( 2 ), 1.0E-3 / equatorialRadius );

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
