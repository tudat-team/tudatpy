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
 *      100903    K. Kumar          File created.
 *      100916    L. Abdulkadir     File checked.
 *      100929    K. Kumar          Checked code by D. Dirkx added.
 *      101110    K. Kumar          Added raiseToExponentPower( ) function.
 *      102410    D. Dirkx          Minor comment changes as code check.
 *      101213    K. Kumar          Bugfix raiseToIntegerExponent( ); renamed raiseToIntegerPower( ).
 *                                  Added computeAbsoluteValue( ) functions.
 *      110202    K. Kumar          Added overload for State* for computeLinearInterpolation( ).
 *      110111    J. Melman         Added computeModulo( ) function.
 *      110411    K. Kumar          Added convertCartesianToSpherical( ) function.
 *      110606    J. Melman         Removed possible singularity from
 *                                  convertCartesianToSpherical.
 *      110707    K. Kumar          Added computeSampleMean( ), computeSampleVariance( ) functions.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120202    K. Kumar          Moved linear interpolation functions into new Interpolators
 *                                  sub-directory.
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of
 *          Scientific Computing. Cambridge University Press, February 2002.
 *
 */

#include "Tudat/Mathematics/BasicMathematics/nearestNeighbourSearch.h"
#include "Tudat/Mathematics/Interpolators/linearInterpolation.h"

namespace tudat
{
namespace mathematics
{
namespace interpolators
{

//! Compute linear interpolation.
double computeLinearInterpolation( Eigen::VectorXd& sortedIndependentVariables,
                                   Eigen::VectorXd& associatedDependentVariables,
                                   double targetIndependentVariableValue )
{
    // Declare local variables.
    // Declare nearest neighbor.
    int nearestNeighbor;
    double locationTargetIndependentVariableValueInInterval;

    // Compute nearest neighbor in sorted vector of independent variables.
    // Result is always to the left of the target independent variable value.
    nearestNeighbor = mathematics
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
Eigen::VectorXd computeLinearInterpolation(
        std::map < double, Eigen::VectorXd >& sortedIndepedentAndDependentVariables,
        double targetIndependentVariableValue )
{
    // Declare local variables.
    // Declare nearest neighbor.
    int nearestLeftNeighbor;

    // Declare location of target independent variable value in interval.
    double locationTargetIndependentVariableValueInInterval;

    // Declare map iterators
    std::map < double, Eigen::VectorXd >::iterator mapIteratorIntervalLeft;
    std::map < double, Eigen::VectorXd >::iterator mapIteratorIntervalRight;

    // Compute nearest neighbor in map of data.
    // Result is always to the left of the target independent variable value.
    nearestLeftNeighbor = mathematics::
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

} // namespace interpolators
} // namespace mathematics
} // namespace tudat
