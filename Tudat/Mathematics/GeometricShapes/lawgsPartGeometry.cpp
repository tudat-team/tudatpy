/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Craidon, C.B. A Desription of the Langley Wireframe Geometry Standard (LaWGS) format, NASA
 *          TECHNICAL MEMORANDUM 85767.
 *
 */

#include <boost/shared_ptr.hpp>
#include "Tudat/Mathematics/GeometricShapes/lawgsPartGeometry.h"

namespace tudat
{
namespace geometric_shapes
{

//! Constructor from surface geometry (i.e., geometry type conversion).
void LawgsPartGeometry::setMesh( boost::shared_ptr< SingleSurfaceGeometry > originalSurface,
                                 int numberOfLinesIn, int numberOfPointsIn )
{
    // Set (temporary name) of part.
    name_ = "copied surface";

    // Set size of mesh.
    numberOfLines_ = numberOfLinesIn;
    numberOfPoints_ = numberOfPointsIn;

    // Allocate mesh points.
    meshPoints_.resize( boost::extents[ numberOfLines_ ][ numberOfPoints_ ] );

    // Set grid sizes from requested number of sample points.
    double independentVariableGridSize1 =
            ( originalSurface->getMaximumIndependentVariable( 1 ) -
              originalSurface->getMinimumIndependentVariable( 1 ) ) /
            ( static_cast< double >( numberOfLines_ - 1 ) );

    double independentVariableGridSize2 =
            ( originalSurface->getMaximumIndependentVariable( 2 ) -
              originalSurface->getMinimumIndependentVariable( 2 ) ) /
            ( static_cast< double >( numberOfPoints_ - 1 ) );

    // Declare sample point variables.
    double variable1, variable2;

    // Declare and set minimum values of independent variables.
    double minimumIndependentVariable1 = originalSurface->getMinimumIndependentVariable( 1 );
    double minimumIndependentVariable2 = originalSurface->getMinimumIndependentVariable( 2 );

    // Loop through the number of lines and points specified and sample
    // geometry at fixed intervals.
    for ( int i = 0; i < numberOfLines_; i++ )
    {
        for ( int j = 0; j < numberOfPoints_; j++ )
        {
            // Set sampling point of original geometry.
            variable1 = minimumIndependentVariable1 + i * independentVariableGridSize1;
            variable2 = minimumIndependentVariable2 + j * independentVariableGridSize2;

            // Set new mesh point.
            meshPoints_[ i ][ j ] = originalSurface->getSurfacePoint( variable1, variable2 );
        }
    }

    // Perform panel calculations for mesh.
    performPanelCalculations( );
}

//! Copy constructor.
LawgsPartGeometry::LawgsPartGeometry( const LawgsPartGeometry& partToCopy )
    : QuadrilateralMeshedSurfaceGeometry( )
{
    // Copy all properties of partToCopy to new part.
    name_ = partToCopy.name_;
    reversalOperator_ = partToCopy.reversalOperator_;
    numberOfLines_ = partToCopy.numberOfLines_;
    numberOfPoints_ = partToCopy.numberOfPoints_;
    rotationMatrix_= partToCopy.rotationMatrix_;
    scalingMatrix_ = partToCopy.scalingMatrix_;
    offset_= partToCopy.offset_;

    // Copy surface points array.
    meshPoints_ = partToCopy.meshPoints_;
}

//! Get surface point.
Eigen::VectorXd LawgsPartGeometry::getSurfacePoint( const double independentVariable1,
                                                    const double independentVariable2 )
{
    // Declare local variables denoting 'start' of panel.
    int pointIndex, lineIndex;

    // Declare local variable denoting surface point.
    Eigen::VectorXd point = Eigen::VectorXd( 3 );

    // Set local variables.
    pointIndex = static_cast< int > ( floor( independentVariable2 ) );
    lineIndex = static_cast< int > ( floor( independentVariable1 ) );

    // Start with panel centroid.
    point = panelCentroids_[ lineIndex ][ pointIndex ];

    // Move back to panel corner.
    point -= 0.5 *( meshPoints_[ lineIndex + 1 ][ pointIndex + 1 ] -
                    meshPoints_[ lineIndex ][ pointIndex ]);

    // Add contribution of 1st independent variable on panel.
    point += ( independentVariable1 - lineIndex ) *
             ( meshPoints_[ lineIndex + 1 ][ pointIndex ]
               - meshPoints_[ lineIndex ][ pointIndex ] );

    // Add contribution of 2nd independent variable on panel.
    point += ( independentVariable2 - pointIndex ) *
             ( meshPoints_[ lineIndex ][ pointIndex + 1 ] -
               meshPoints_[ lineIndex ][ pointIndex ] );

    return point;
}

//! Get surface derivative (currently not implemented).
Eigen::VectorXd LawgsPartGeometry::getSurfaceDerivative( const double u, const double v,
                                                         const int uDerivative,
                                                         const int vDerivative )
{
    std::string errorMessage =  "Warning, surface derivative function not implemented in LawgsPartGeometry class. Not able to return the "
            + std::to_string( uDerivative ) + ", "
            + std::to_string( vDerivative ) + "the derivative at point,"
            + std::to_string( u ) + ", " + std::to_string( v );
    throw std::runtime_error( errorMessage );

    return Eigen::Vector3d( 0.0, 0.0, 0.0 );
}

//! Get parameter.
double LawgsPartGeometry::getParameter( const int parameterIndex )
{
    std::string errorMessage =  "WWarning, get parameter function not implemented in LawgsPartGeometry class, unable to retrieve parameter"
            + std::to_string( parameterIndex );
    throw std::runtime_error( errorMessage );

    return 0.0;
}

//! Set parameter.
void LawgsPartGeometry::setParameter( const int parameterIndex, const double value )
{
    std::string errorMessage =  "WWarning, set parameter function not implemented in LawgsPartGeometry class, unable to set parameter"
            + std::to_string( parameterIndex );
    throw std::runtime_error( errorMessage );
}

//! Overload ostream to print class information.
std::ostream& operator << ( std::ostream& stream, LawgsPartGeometry& lawgsPartGeometry )
{
    stream << "This is a Langley Wireframe Geometry Standard surface geometry"
           << " of a single part." << std::endl;
    stream << "The number of lines ( contours ) is: "
           << lawgsPartGeometry.numberOfLines_ << std::endl;
    stream << "The number of points per line is: "
           << lawgsPartGeometry.numberOfPoints_ << std::endl;
    stream << "The part name is: " << lawgsPartGeometry.name_ << std::endl;

    // Return stream.
    return stream;
}

} // namespace geometric_shapes
} // namespace tudat
