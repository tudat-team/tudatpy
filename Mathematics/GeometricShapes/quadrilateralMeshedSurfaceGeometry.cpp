/*! \file quadrilateralMeshedSurfaceGeometry.cpp
 *  This file contains the definition of the Quadrilateral Meshed Surface
 *  Geometry class
 *
 *  Path              : Mathematics/Geometry/
 *  Version           : 1
 *  Check status      : Unchecked
 *
 *  Author            : Dominic Dirkx
 *  Affiliation       : TU Delft
 *  E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *  Checker           : J. Melman
 *  Affiliation       : Delft University of Technology
 *  E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *  Date created      : 25 November, 2010
 *  Last modified     : 6 February, 2011
 *
 *  References
 *    An example of a heritage code which uses such a mesh is found in:
 *    The Mark IV Supersonic-Hypersonic Arbitrary Body Program, Volume II-
 *    Program Formulation, Douglas Aircraft Company, AFFDL-TR-73-159, Volume II.
 *
 *  Notes
 *
 *    This class uses pointers to pointers to denote two-dimensional arrays
 *    for heritage reasons, instead of the more modern C++ vector or map types.
 *
 *    The numberOfLines_ and numberOfPoints_ member variables denote the number of
 *    mesh points. The number of panels in the mesh will be numberOfLines_  - 1 by
 *    numberOfPoints_ -1.
 *
 *  Copyright (c) 2010 Delft University of Technology.
 *
 *  This software is protected by national and international copyright.
 *  Any unauthorized use, reproduction or modification is unlawful and
 *  will be prosecuted. Commercial and non-private application of the
 *  software in any form is strictly prohibited unless otherwise granted
 *  by the authors.
 *
 *  The code is provided without any warranty; without even the implied
 *  warranty of merchantibility or fitness for a particular purpose.
 *
 *  Changelog
 *    YYMMDD   Author               Comment
 *    101125   D. Dirkx             First version of file
 *    110127   D. Dirkx             Finalized for code check.
 *    110206   J. Melman            Minor formatting issues.
 */

// Include statements.
#include "quadrilateralMeshedSurfaceGeometry.h"

//! Default constructor
QuadrilateralMeshedSurfaceGeometry::QuadrilateralMeshedSurfaceGeometry( ) :
        numberOfLines_( -0 ), numberOfPoints_( -0 ), reversalOperator_( 1 ),
        meshPoints_( NULL ), panelCentroids_( NULL ),
        panelSurfaceNormals_( NULL ), panelAreas_( NULL ), totalArea_( -0.0 )
{
}

//! Default destructor
QuadrilateralMeshedSurfaceGeometry::~QuadrilateralMeshedSurfaceGeometry( )
{
    // Loop through all lines on surface to deallocate variables. Set pointers
    // to NULL afterwards to prevent possible problems when calling this
    // function multiple times.
    for ( int j = 0; j < numberOfLines_ - 1; j++)
    {
        delete[ ] meshPoints_[ j ];
        delete[ ] panelCentroids_[ j ];
        delete[ ] panelSurfaceNormals_[ j ];
        delete[ ] panelAreas_[ j ];
        meshPoints_[ j ] = NULL;
        panelCentroids_[ j ] = NULL;
        panelSurfaceNormals_[ j ] = NULL;
        panelAreas_[ j ] = NULL;
    }

    // Delete final surface point line.
    delete[ ] meshPoints_[ numberOfLines_ - 1 ];
    meshPoints_[ numberOfLines_ - 1 ] = NULL;

    // Delete entire pointers and reset to NULL.
    delete[ ] meshPoints_;
    delete[ ] panelCentroids_;
    delete[ ] panelSurfaceNormals_;
    delete[ ] panelAreas_;
    meshPoints_ = NULL;
    panelCentroids_ = NULL;
    panelSurfaceNormals_ = NULL;
    panelAreas_ = NULL;
}

//! Function to calculate panel characteristics
void QuadrilateralMeshedSurfaceGeometry::performPanelCalculations( )
{
    // Allocate memory for panel properties.
    panelCentroids_ = new Vector3d* [ numberOfLines_ - 1 ];
    panelSurfaceNormals_ = new Vector3d* [ numberOfLines_ - 1 ];
    panelAreas_ = new double* [ numberOfLines_ - 1 ];

    // Allocate memory for panel properties per line.
    for ( int i = 0; i < numberOfLines_ - 1 ; i++ )
    {
        panelCentroids_[ i ] = new Vector3d[ numberOfPoints_ - 1 ];
        panelSurfaceNormals_[ i ] = new Vector3d[ numberOfPoints_ - 1 ];
        panelAreas_[ i ] = new double[ numberOfPoints_ - 1 ];
    }

    // Declare local variables for normal and area determination.
    Vector3d crossVector1;
    Vector3d crossVector2;

    // Reset total area.
    totalArea_ = 0.0;

    // Loop over all panels to determine properties.
    for ( int i = 0; i < numberOfLines_ - 1; i++ )
    {
        for ( int j = 0; j < numberOfPoints_ - 1; j++ )
        {
            // Set panel centroid.
            panelCentroids_[ i ][ j ] = ( meshPoints_[ i ][ j ] +
                                          meshPoints_[ i + 1 ][ j ] +
                                          meshPoints_[ i ][ j + 1 ] +
                                          meshPoints_[ i + 1 ][ j + 1 ] ) / 4;

            // Set panel cross vectors.
            crossVector1 = meshPoints_[ i + 1 ][ j + 1 ] - meshPoints_[ i ][ j ];
            crossVector2 = meshPoints_[ i + 1 ][ j ] - meshPoints_[ i ][ j + 1 ];

            // Set panel normal (not yet normalized).
            panelSurfaceNormals_[ i ][ j ] = crossVector1.cross( crossVector2 );

            // Set panel area (not yet correct size).
            panelAreas_[ i ][ j ] = panelSurfaceNormals_[ i ][ j ].norm( );
            if ( panelAreas_[ i ][ j ] < mathematics::MACHINE_PRECISION_DOUBLES )
            {
                std::cerr << "WARNING panel area is zero in part at panel" << i
                          << ", " << j << std::endl;               
            }            

            // Normalize panel normal and, if necessary, invert normal direction.
            panelSurfaceNormals_[ i ][ j ] *= reversalOperator_;
            panelSurfaceNormals_[ i ][ j ].normalize( );

            // Set panel area to correct size.
            panelAreas_[ i ][ j ] *= 0.5;

            // Add panel area to total area.
            totalArea_ += panelAreas_[ i ][ j ];
        }
    }
}

//! Function to retrieve a surface point.
Vector3d QuadrilateralMeshedSurfaceGeometry::getMeshPoint(
        const int& lineIndex,
        const int& pointIndex )
{
    return meshPoints_[ lineIndex ][ pointIndex ];
}

//! Function to retrieve a panel area.
double QuadrilateralMeshedSurfaceGeometry::getPanelArea(
        const int& lineIndex,
        const int& pointIndex )
{
    return panelAreas_[ lineIndex ][ pointIndex ];
}

//! Function to retrieve a panel centroid.
Vector3d QuadrilateralMeshedSurfaceGeometry::getPanelCentroid(
        const int& lineIndex,
        const int& pointIndex )
{
    return panelCentroids_[ lineIndex ][ pointIndex ];
}

//! Function to retrieve an outward panel surface normal
Vector3d QuadrilateralMeshedSurfaceGeometry::getPanelSurfaceNormal(
        const int& lineIndex,
        const int& pointIndex )
{
    return panelSurfaceNormals_[ lineIndex ][ pointIndex ];
}

//! Function to retrieve number of lines.
int QuadrilateralMeshedSurfaceGeometry::getNumberOfLines( )
{
    return numberOfLines_;
}

//! Function to retrieve number of points.
int QuadrilateralMeshedSurfaceGeometry::getNumberOfPoints( )
{
    return numberOfPoints_;
}

//! Returns the total area of the mesh.
double QuadrilateralMeshedSurfaceGeometry::getTotalArea( )
{
    return totalArea_;
}

//! Sets reversal operator.
void QuadrilateralMeshedSurfaceGeometry::setReversalOperator(
        bool isMeshInverted )
{
    if( isMeshInverted == 0 )
    {
        reversalOperator_ = 1;
    }
    else
    {
        reversalOperator_ = 0;
    }
}

//! Returns a boolean denoting if the mesh is inverted.
 bool QuadrilateralMeshedSurfaceGeometry::getReversalOperator( )
 {
     bool isMeshInverted;
     if( reversalOperator_ == 1 )
     {
         isMeshInverted = 0;
     }
     else
     {
         isMeshInverted = 1;
     }
     return isMeshInverted;
 }

 //! Overload ostream to print class information.
 std::ostream& operator<<( std::ostream& stream,
                           QuadrilateralMeshedSurfaceGeometry& geometry )
 {
     stream << "This is a quadrilateral meshed surface geometry"
            << " of a single part." << std::endl;
     stream << "The number of lines ( contours ) is: "
            << geometry.numberOfLines_ << std::endl;
     stream << "The number of points per line is: "
            << geometry.numberOfPoints_ << std::endl;

     // Return stream.
     return stream;
 }

 // End of file.
