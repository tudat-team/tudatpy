/*! \file basicMathematicsFunctions.cpp
 *    Source file that defines the basicMathematicsFunctions namespace,
 *    containing all basic functions contained in Tudat.
 *
 *    Path              : /Mathematics/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Author            : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Checker           : L. Abdulkadir
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : L.Abdulkadir@student.tudelft.nl
 *
 *    Date created      : 3 september, 2010
 *    Last modified     : 29 september, 2010
 *
 *    References
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
 *      YYMMDD    author              comment
 *      100903    K. Kumar            File created.
 *      100916    L. Abdulkadir       File checked.
 *      100929    K. Kumar            Checked code by D. Dirkx added.
 */

// Include statements.
#include "basicMathematicsFunctions.h"

//! mathematics namespace.
namespace mathematics
{
//! Linear interpolation.
double computeLinearInterpolation( Vector& sortedIndependentVariables,
                                   Vector& associatedDependentVariables,
                                   double& targetIndependentVariableValue )
{

    // Declare local variables.
    // Declare nearest neighbor.
    int nearestNeighbor;
    double locationTargetIndependentVariableValueInInterval;

    // Compute nearest neighbor in sorted vector of independent variables.
    // Result is always to the left of the target independent variable value.
    nearestNeighbor = basic_functions
                      ::computeNearestLeftNeighborUsingBinarySearch(
            sortedIndependentVariables,
            targetIndependentVariableValue );

    // Compute location of target independent variable value in interval
    // between nearest neighbors.
    locationTargetIndependentVariableValueInInterval
            = ( targetIndependentVariableValue
              - sortedIndependentVariables[ nearestNeighbor ] )
             / ( sortedIndependentVariables[ nearestNeighbor + 1 ]
                 - sortedIndependentVariables[ nearestNeighbor ] );

    // Return the computed value of the dependent variable.
    return ( associatedDependentVariables[ nearestNeighbor ]
             * ( 1 - locationTargetIndependentVariableValueInInterval )
             + associatedDependentVariables[ nearestNeighbor + 1 ]
             * locationTargetIndependentVariableValueInInterval );
}

//! Linear interpolation.
Vector computeLinearInterpolation(
        std::map < double, Vector >& sortedIndepedentAndDependentVariables,
        double& targetIndependentVariableValue )
{
    // Declare local variables.
    // Declare nearest neighbor.
    int nearestLeftNeighbor;

    // Declare location of target independent variable value in interval.
    double locationTargetIndependentVariableValueInInterval;

    // Declare map iterators
    std::map < double, Vector >::iterator mapIteratorIntervalLeft;
    std::map < double, Vector >::iterator mapIteratorIntervalRight;

    // Compute nearest neighbor in map of data.
    // Result is always to the left of the target independent variable value.
    nearestLeftNeighbor = basic_functions::
                          computeNearestLeftNeighborUsingBinarySearch(
                                  sortedIndepedentAndDependentVariables,
                                  targetIndependentVariableValue );

    // Compute location of target independent variable value in interval
    // between nearest neighbors.
    mapIteratorIntervalLeft = sortedIndepedentAndDependentVariables.begin( );
    advance( mapIteratorIntervalLeft, nearestLeftNeighbor );
    mapIteratorIntervalRight = sortedIndepedentAndDependentVariables.begin( );
    advance( mapIteratorIntervalRight, nearestLeftNeighbor + 1 );
    locationTargetIndependentVariableValueInInterval
            = ( targetIndependentVariableValue
              - mapIteratorIntervalLeft->first )
             / ( mapIteratorIntervalRight->first
                 - mapIteratorIntervalLeft->first );

    // Return the computed value of the dependent variable.
    return ( mapIteratorIntervalLeft->second
             * ( 1 - locationTargetIndependentVariableValueInInterval )
             + mapIteratorIntervalRight->second
             * locationTargetIndependentVariableValueInInterval );
}

//! Converts spherical to cartesian coordinates.
void convertSphericalToCartesian( const double& radius,
                                  const double& azimuthAngle,
                                  const double& zenithAngle,
                                  Vector& cartesianCoordinates )
{
    // Declaring sine and cosine which have multiple usages to save computation
    // time.
    double cosineOfAzimuthAngle = cos ( azimuthAngle );
    double sineOfZenithAngle = sin( zenithAngle );

    // Perform transformation.
    cartesianCoordinates( 0 ) = radius * cosineOfAzimuthAngle
                                * sineOfZenithAngle;
    cartesianCoordinates( 1 ) = radius * sin( azimuthAngle )
                                * sineOfZenithAngle;
    cartesianCoordinates( 2 ) = radius * cos( zenithAngle );
}

//! Converts cylindrical to cartesian coordinates, z value left unaffected.
void convertCylindricalToCartesian( const double& radius,
                                    const double& azimuthAngle,
                                    Vector& cartesianCoordinates )
{
    // Perform transformation, z value should be set outside function.
    cartesianCoordinates( 0 ) = radius * cos( azimuthAngle );
    cartesianCoordinates( 1 ) = radius * sin( azimuthAngle );
}

}

// End of file.
