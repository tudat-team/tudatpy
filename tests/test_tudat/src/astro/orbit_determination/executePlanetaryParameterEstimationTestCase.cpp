/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "tudat/simulation/estimation_setup/executePlanetaryParameterEstimationTestCase.h"

namespace tudat
{
namespace unit_tests
{

#if TUDAT_BUILD_EXPLICIT_INSTANTIATIONS
template std::pair< std::shared_ptr< simulation_setup::EstimationOutput< double, double > >, Eigen::VectorXd >
executePlanetaryParameterEstimation< double, double >( const int observableType,
                                                       Eigen::VectorXd parameterPerturbation,
                                                       Eigen::MatrixXd inverseAPrioriCovariance,
                                                       const double weight );
#endif

}  // namespace unit_tests

}  // namespace tudat
