/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PROPAGATECOVARIANCE_H
#define TUDAT_PROPAGATECOVARIANCE_H

#include<Tudat/Astrodynamics/Propagators/stateTransitionMatrixInterface.h>

namespace tudat
{

namespace propagators
{


std::map< double, Eigen::MatrixXd > propagateCovariance(
        const Eigen::MatrixXd& initialCovariance,
        const std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const double timeStep,
        const double initialTime,
        const double finalTime );

std::map< double, Eigen::VectorXd > propagateFormalErrors(
        const Eigen::MatrixXd& initialCovariance,
        const std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const double timeStep,
        const double initialTime,
        const double finalTime );

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATECOVARIANCE_H
