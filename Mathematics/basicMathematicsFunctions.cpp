/*! \file basicMathematicsFunctions.cpp
 *    Source file that defines the basicMathematicsFunctions namespace,
 *    containing all basic functions contained in Tudat.
 *
 *    Path              : /Mathematics/
 *    Version           : 12
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Author            : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Author            : J.C.P Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Checker           : L. Abdulkadir
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : L.Abdulkadir@student.tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Date created      : 3 September, 2010
 *    Last modified     : 5 September, 2011
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of
 *          Scientific Computing. Cambridge University Press, February 2002.
 *
 *    Notes
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      100903    K. Kumar          File created.
 *      100916    L. Abdulkadir     File checked.
 *      100929    K. Kumar          Checked code by D. Dirkx added.
 *      101110    K. Kumar          Added raiseToExponentPower() function.
 *      102410    D. Dirkx          Minor comment changes as code check.
 *      101213    K. Kumar          Bugfix raiseToIntegerExponent(); renamed raiseToIntegerPower().
 *                                  Added computeAbsoluteValue() functions.
 *      110202    K. Kumar          Added overload for State* for computeLinearInterpolation().
 *      110111    J. Melman         Added computeModulo() function.
 *      110411    K. Kumar          Added convertCartesianToSpherical() function.
 *      110606    J. Melman         Removed possible singularity from
 *                                  convertCartesianToSpherical.
 *      110707    K. Kumar          Added computeSampleMean(), computeSampleVariance() functions.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 */

// Include statements.
#include <cmath>
#include <iostream>
#include <numeric>
#include "Basics/basicFunctions.h"
#include "Mathematics/basicMathematicsFunctions.h"

// Using declarations.
using std::map;
using std::vector;
using std::accumulate;
using std::pow;
using mathematics::MACHINE_PRECISION_DOUBLES;

//! Mathematics namespace.
namespace mathematics
{

//! Compute linear interpolation.
double computeLinearInterpolation( VectorXd& sortedIndependentVariables,
                                   VectorXd& associatedDependentVariables,
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

//! Compute linear interpolation.
VectorXd computeLinearInterpolation(
        std::map < double, VectorXd >& sortedIndepedentAndDependentVariables,
        double& targetIndependentVariableValue )
{
    // Declare local variables.
    // Declare nearest neighbor.
    int nearestLeftNeighbor;

    // Declare location of target independent variable value in interval.
    double locationTargetIndependentVariableValueInInterval;

    // Declare map iterators
    std::map < double, VectorXd >::iterator mapIteratorIntervalLeft;
    std::map < double, VectorXd >::iterator mapIteratorIntervalRight;

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

//! Compute linear interpolation.
State* computeLinearInterpolation(
        std::map < double, State* >& sortedIndepedentAndDependentVariables,
        double& targetIndependentVariableValue )
{
    // Declare local variables.
    // Declare nearest neighbor.
    int nearestLeftNeighbor;

    // Pointer to state.
    State* pointerToState_ = new State;

    // Declare location of target independent variable value in interval.
    double locationTargetIndependentVariableValueInInterval;

    // Declare map iterators
    std::map < double, State* >::iterator mapIteratorIntervalLeft;
    std::map < double, State* >::iterator mapIteratorIntervalRight;

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

    // Set vector in pointer to State object to the computed value of the
    // dependent variable.

    pointerToState_->state
            = mapIteratorIntervalLeft->second->state
              * ( 1 - locationTargetIndependentVariableValueInInterval )
              + mapIteratorIntervalRight->second->state
              * locationTargetIndependentVariableValueInInterval;

    // Return pointer to State object.
    return pointerToState_;
}

//! Convert spherical to cartesian coordinates.
void convertSphericalToCartesian( const double& radius,
                                  const double& azimuthAngle,
                                  const double& zenithAngle,
                                  VectorXd& cartesianCoordinates )
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

//! Convert cartesian to spherical coordinates.
void convertCartesianToSpherical( const VectorXd& cartesianCoordinates,
                                  VectorXd& sphericalCoordinates )
{
    // Compute transformation of Cartesian coordinates to spherical
    // coordinates.
    sphericalCoordinates( 0 ) = cartesianCoordinates.norm( );

    // Check if coordinates are at origin.
    if ( sphericalCoordinates( 0 ) < MACHINE_PRECISION_DOUBLES )
    {
        sphericalCoordinates( 1 ) = 0.0;
        sphericalCoordinates( 2 ) = 0.0;
    }

    // Else compute coordinates using trigonometric relationships.
    else
    {
        sphericalCoordinates( 1 ) = atan2( cartesianCoordinates( 1 ),
                                       cartesianCoordinates( 0 ) );
        sphericalCoordinates( 2 ) = acos( cartesianCoordinates( 2 )
                                      / sphericalCoordinates( 0 ) );
    }
}

//! Convert cylindrical to cartesian coordinates, z value left unaffected.
void convertCylindricalToCartesian( const double& radius,
                                    const double& azimuthAngle,
                                    VectorXd& cartesianCoordinates )
{
    // Perform transformation, z value should be set outside function.
    cartesianCoordinates( 0 ) = radius * cos( azimuthAngle );
    cartesianCoordinates( 1 ) = radius * sin( azimuthAngle );
}

//! Compute modulo of double.
double computeModulo( const double& dividend, const double& divisor )
{
    return dividend - divisor * floor( dividend / divisor );
}

//! Compute sample mean.
double computeSampleMean( const vector< double >& sampleData )
{
    // Return sample mean.
    return accumulate( sampleData.begin( ), sampleData.end( ), 0.0 )
            / static_cast< double >( sampleData.size( ) );
}

//! Compute sample variance.
double computeSampleVariance( const vector< double >& sampleData )
{
    // Declare local variables.
    // Declare and compute sample mean.
    double sampleMean_ = computeSampleMean( sampleData );

    // Declare and initialize sum of residuals squared.
    double sumOfResidualsSquared_ = 0.0;

    // Compute sum of residuals of sample data squared.
    for ( unsigned int i = 0; i < sampleData.size( ); i++ )
    {
        sumOfResidualsSquared_ += pow( sampleData.at( i ) - sampleMean_, 2.0 );
    }

    // Return sample variance.
    return 1.0 / ( static_cast< double >( sampleData.size( ) ) - 1.0 )
            * sumOfResidualsSquared_;
}

}

// End of file.
