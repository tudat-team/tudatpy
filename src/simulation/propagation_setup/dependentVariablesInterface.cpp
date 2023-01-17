/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/propagation_setup/dependentVariablesInterface.h"

namespace tudat
{

namespace propagators
{



////! Function to reset the dependent variables interpolators
//template< typename TimeType >
//void MultiArcDependentVariablesInterface< TimeType >::updateDependentVariablesInterpolators(
//        const std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, Eigen::VectorXd > > > dependentVariablesInterpolators,
//        const std::vector< double >& arcStartTimes,
//        const std::vector< double >& arcEndTimes )
//{
//
//}

////! Function to get the single-arc dependent variable at a given time.
//template< typename TimeType >
//Eigen::VectorXd MultiArcDependentVariablesInterface< TimeType >::getDependentVariables(
//        const TimeType evaluationTime )
//{
//
//}


////! Function to retrieve the current arc for a given time
//template< typename TimeType >
//std::pair< int, double > MultiArcDependentVariablesInterface< TimeType >::getCurrentArc( const TimeType evaluationTime )
//{
//
//
//}


////! Function to get the dependent variable at a given time.
//template< typename TimeType >
//Eigen::VectorXd HybridArcDependentVariablesInterface< TimeType >::getDependentVariables(
//        const TimeType evaluationTime )
//{
//
//}


} // namespace propagators

} // namespace tudat
