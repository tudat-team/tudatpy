/*    Copyright (c) 2010-2013, Delft University of Technology
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
 *      110206    D. Dirkx          First version of file.
 *      110208    D. Dirkx          Finalized for code check.
 *      110212    J. Melman         Fixed many minor formatting issues.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      120628    A. Ronse          Boostified unit test.
 *
 *    References
 *      Craidon, C.B. A Desription of the Langley Wireframe Geometry Standard (LaWGS) format, NASA
 *          TECHNICAL MEMORANDUM 85767.
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <Eigen/Core>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Mathematics/GeometricShapes/lawgsPartGeometry.h"
#include "Tudat/Mathematics/GeometricShapes/sphereSegment.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_Lawgs_Surface_Geometry )

//! Test implementation of Lawgs surface geometry.
BOOST_AUTO_TEST_CASE( testLawgsSurfaceGeometry )
{
    using namespace tudat;
    using namespace tudat::geometric_shapes;

    // Create a full sphere as test geometry, with a radius of 2.0.
    const double sphereRadius = 2.0;
    boost::shared_ptr< SphereSegment > sphere = boost::make_shared< SphereSegment >(
                sphereRadius );

    // Create a Lawgs mesh of the sphere.
    LawgsPartGeometry lawgsSurface;
    const int numberOfLines = 21;
    const int numberOfPoints = 21;
    lawgsSurface.setMesh( sphere, numberOfLines, numberOfPoints );

    // Retrieve the total surface area and check if it is sufficiently close
    // to the expected value.
    using tudat::basic_mathematics::mathematical_constants::PI;
    const double totalArea = lawgsSurface.getTotalArea( );
    BOOST_CHECK_SMALL( std::fabs( totalArea - 4.0 * PI
                                  * ( std::pow( sphereRadius, 2.0 ) ) ), 0.6 );

    // Test if number of lines on mesh is correct.
    BOOST_CHECK_EQUAL( lawgsSurface.getNumberOfLines( ), numberOfLines );

    // Test if number of points per line on mesh is correct.
    BOOST_CHECK_EQUAL( lawgsSurface.getNumberOfPoints( ), numberOfPoints );

    // Set part name.
    std::string partName = "sphere";
    lawgsSurface.setName( partName );

    // Test if part name is properly retrieved.
    BOOST_CHECK_EQUAL( lawgsSurface.getName( ), partName );

    // Retrieve normal and centroid for panel 0, 0.
    Eigen::Vector3d testNormal = lawgsSurface.getPanelSurfaceNormal( 0, 0 );
    Eigen::Vector3d testCentroid = lawgsSurface.getPanelCentroid( 0, 0 );

    // Test whether centroid and normal are collinear for panel 0, 0.
    BOOST_CHECK_SMALL( std::fabs( testCentroid.normalized( ).dot(
                                      testNormal.normalized( ) ) ) - 1.0, 1.0e-5 );

    // Test if the position of the x- and y-coordinate of panel 0, 0 is correct.
    BOOST_CHECK_SMALL( std::fabs( std::atan( testCentroid.y( ) / testCentroid.x( ) ) - PI / 20.0 ),
                       std::numeric_limits< double >::epsilon( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
