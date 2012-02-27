/*    Copyright (c) 2010-2012 Delft University of Technology.
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
 *      YYMMDD    Author            Comment
 *      101125    D. Dirkx          First version of file.
 *      110127    D. Dirkx          Finalized for code check.
 *      110206    J. Melman         Minor formatting issues. Identified multiple warnings.
 *      110207    D. Dirkx          Fixed warning problems by extending cerr comments.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *
 *    References
 *      Craidon, C.B. A Desription of the Langley Wireframe Geometry Standard
 *          (LaWGS) format, NASA TECHNICAL MEMORANDUM 85767.
 *
 */

#include "Tudat/Mathematics/GeometricShapes/lawgsPartGeometry.h"

namespace tudat
{

//! Constructor from surface geometry (i.e., geometry type conversion).
void LawgsPartGeometry::setMesh( SingleSurfaceGeometry* originalSurface,
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
Eigen::VectorXd LawgsPartGeometry::getSurfacePoint( double independentVariable1,
                                                    double independentVariable2 )
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
Eigen::VectorXd LawgsPartGeometry::getSurfaceDerivative( double u, double v,
                                                         int uDerivative, int vDerivative )
{
    std::cerr << "Surface derivative function not implemented in "
              << "LawgsPartGeometry class. Not able to return the "
              << uDerivative << ", " << vDerivative << "th derivative at point,"
              << u << ", " << v << ". Returning zero vector." << std::endl;

    return Eigen::Vector3d( 0.0, 0.0, 0.0 );
}

//! Get parameter.
double LawgsPartGeometry::getParameter( int parameterIndex )
{
    std::cerr << "Get parameter function not implemented in LawgsPartGeometry"
              << "class, unable to retrieve parameter "<< parameterIndex
              << ". Returning zero." << std::endl;

    return 0.0;
}

//! Set parameter.
void LawgsPartGeometry::setParameter( int parameterIndex, double value )
{
    std::cerr << "Set parameter function not implemented in LawgsPartGeometry"
              << "class. Unable to set value of " << value << " at parameter index "
              << parameterIndex << std::endl;
}

//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream,
                          LawgsPartGeometry& lawgsPartGeometry )
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

} // namespace tudat
