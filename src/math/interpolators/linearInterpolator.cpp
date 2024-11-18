/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of Scientific Computing. Cambridge
 *          University Press, February 2002.
 *
 */

#include <boost/multi_array.hpp>

#include "tudat/math/interpolators/linearInterpolator.h"

namespace tudat
{
namespace interpolators
{

//! Compute linear interpolation.
double computeLinearInterpolation( const Eigen::VectorXd& sortedIndependentVariables,
                                   const Eigen::VectorXd& associatedDependentVariables,
                                   const double targetIndependentVariableValue )
{
    // Declare local variables.
    // Declare nearest neighbor.
    int nearestNeighbor;
    double locationTargetIndependentVariableValueInInterval;

    // Compute nearest neighbor in sorted vector of independent variables.
    // Result is always to the left of the target independent variable value.
    nearestNeighbor = basic_mathematics::computeNearestLeftNeighborUsingBinarySearch(
            sortedIndependentVariables, targetIndependentVariableValue );

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
        const std::map < double, Eigen::VectorXd >& sortedIndepedentAndDependentVariables,
        const double targetIndependentVariableValue )
{
    // Declare local variables.
    // Declare nearest neighbor.
    int nearestLeftNeighbor;

    // Declare location of target independent variable value in interval.
    double locationTargetIndependentVariableValueInInterval;

    // Declare map iterators
    std::map< double, Eigen::VectorXd >::const_iterator mapIteratorIntervalLeft;
    std::map< double, Eigen::VectorXd >::const_iterator mapIteratorIntervalRight;

    // Compute nearest neighbor in map of data.
    // Result is always to the left of the target independent variable value.
    nearestLeftNeighbor = basic_mathematics::computeNearestLeftNeighborUsingBinarySearch(
                sortedIndepedentAndDependentVariables, targetIndependentVariableValue );

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

//template class LinearInterpolator< double, Eigen::VectorXd >;
//template class LinearInterpolator< double, Eigen::Vector6d >;
//template class LinearInterpolator< double, Eigen::MatrixXd >;

//template class LinearInterpolator< double, Eigen::Matrix< long double, Eigen::Dynamic, 1 > >;
//template class LinearInterpolator< double, Eigen::Matrix< long double, Eigen::Dynamic, 6 > >;
//template class LinearInterpolator< double, Eigen::Matrix< long double, Eigen::Dynamic,  Eigen::Dynamic > >;

} // namespace interpolators
} // mamespace tudat
