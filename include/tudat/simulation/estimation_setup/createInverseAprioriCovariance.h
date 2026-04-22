/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATE_INVERSE_APRIORI_COVARIANCE_H
#define TUDAT_CREATE_INVERSE_APRIORI_COVARIANCE_H

#include <memory>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameterSet.h"

namespace tudat
{
namespace simulation_setup
{
std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > aprioriUncertaintyEntry(
        const estimatable_parameters::EstimatebleParameterIdentifier& parameterIdentifier,
        const double uncertainty );

std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > aprioriUncertaintyEntry(
        const estimatable_parameters::EstimatebleParameterIdentifier& parameterIdentifier,
        const Eigen::VectorXd& uncertaintyValues );

Eigen::MatrixXd addInverseAprioriCovarianceEntries(
        const Eigen::MatrixXd& inverseAprioriCovariance,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > >& parameterSet,
        const std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > >&
                aprioriUncertaintyPerParameter );

Eigen::MatrixXd createInverseAprioriCovariance(
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > >& parameterSet,
        const std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > >&
                aprioriUncertaintyPerParameter );

}  // namespace simulation_setup
}  // namespace tudat

#endif  // TUDAT_CREATE_INVERSE_APRIORI_COVARIANCE_H
