/*! \file unitTestLawgsSurfaceGeometry.cpp
 *    This file contains the definition of the LawgsPartGeometry unit test.
 *
 *    Path              : Mathematics/Geometry/
 *    Version           : 3
 *    Check status      : Unchecked
 *
 *    Author            : Dominic Dirkx
 *    Affiliation       : TU Delft
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 6 February, 2011
 *    Last modified     : 12 February, 2011
 *
 *    References
 *      C.B. Craidon, "A Description of the Langley Wireframe Geometry Standard
 *      (LaWGS) format", NASA TECHNICAL MEMORANDUM 85767.
 *
 *    Notes
 *
 *    Copyright (c) 2010 Delft University of Technology.
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
 *      YYMMDD   Author               Comment
 *      110206   D. Dirkx             First version of file.
 *      110208   D. Dirkx             Finalized for code check.
 *      110212   J. Melman            Fixed many minor formatting issues.
 */

// Include statements.
#include "unitTestLawgsSurfaceGeometry.h"

//! Namespace for all unit tests.
bool unit_tests::testLawgsSurfaceGeometry( )
{
    // Test result initialised to false.
    bool isLawgsSurfaceGeometryBad_ = false;

    // Create a full sphere as test geometry.
    SphereSegment sphere = SphereSegment();
    double sphereRadius_ = 2.0;
    sphere.setRadius( sphereRadius_ );
    sphere.setMaximumAzimuthAngle( 2.0 * M_PI );
    sphere.setMinimumAzimuthAngle( 0.0 );
    sphere.setMaximumZenithAngle( M_PI );
    sphere.setMinimumZenithAngle( 0.0 );

    // Create a Lawgs mesh of the sphere.
    LawgsPartGeometry lawgsSurface;
    int numberOfLines = 21;
    int numberOfPoints = 21;
    lawgsSurface.setMesh( &sphere, numberOfLines, numberOfPoints );

    // Retrieve the total surface area and check if it is sufficiently close
    // to the expected value.
    double totalArea_ = lawgsSurface.getTotalArea( );

    if( mathematics::computeAbsoluteValue( totalArea_ - 4.0 * M_PI * (
            mathematics::raiseToIntegerPower( sphereRadius_, 2 ) ) ) > 0.6 )
    {
        std::cerr << "Total mesh area does not match sphere area "
                  << "sufficiently well." << std::endl;
        isLawgsSurfaceGeometryBad_ = true;
    }

    // Test if number of lines on mesh is correct.
    if( lawgsSurface.getNumberOfLines( ) != numberOfLines )
    {
        std::cerr << " Number of lines in mesh incorrect." << std::endl;
        isLawgsSurfaceGeometryBad_ = true;
    }

    // Test if number of points per line on mesh is correct.
    if( lawgsSurface.getNumberOfPoints( ) != numberOfPoints )
    {
        std::cerr << " Number of points in mesh is incorrect." << std::endl;
        isLawgsSurfaceGeometryBad_ = true;
    }

    // Set part name.
    std::string partName_ = "sphere";
    lawgsSurface.setName( partName_ );

    // Test if part name is properly retrieved.
    if ( lawgsSurface.getName() != partName_ )
    {
        std::cerr << " Error in part name of mesh." << std::endl;
        isLawgsSurfaceGeometryBad_ = true;
    }

    // Retrieve normal and centroid for panel 0, 0.
    Vector3d testNormal_ = lawgsSurface.getPanelSurfaceNormal( 0, 0 );
    Vector3d testCentroid_ = lawgsSurface.getPanelCentroid( 0, 0 );

    // Test whether centroid and normal are collinear for panel 0, 0.
    if( mathematics::computeAbsoluteValue( testCentroid_.normalized( ).dot(
            testNormal_.normalized( ) ) - 1.0 ) > 1.0e-5 )
    {
        std::cerr << "Normal and centroid of sphere segment mesh not collinear."
                  << std::endl;
        isLawgsSurfaceGeometryBad_ = true;
    }

    // Test if the position of the x- and y-coordinate of panel 0, 0 is correct.
    if( mathematics::computeAbsoluteValue( atan( testCentroid_.y( ) /
            testCentroid_.x( ) ) - M_PI / 20.0 ) >
        mathematics::MACHINE_PRECISION_DOUBLES )
    {
        std::cerr << "x- and y-coordinate of centroid of panel 0, 0 of "
                  << "sphere mesh is incorrect." << std::endl;
        isLawgsSurfaceGeometryBad_ = true;
    }

    // Return test result.
    // If test is succesful return 0, if test fails, return 1.
    return  isLawgsSurfaceGeometryBad_;
}

// End of file.
