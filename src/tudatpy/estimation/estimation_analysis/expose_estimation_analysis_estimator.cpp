/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_estimation_analysis_estimator.h"

#include <pybind11/pybind11.h>

#include "expose_estimation_analysis_estimator_bindings.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"

namespace py = pybind11;

namespace tudatpy
{
namespace estimation
{
namespace estimation_analysis
{

void expose_estimation_analysis_estimator( py::module& m )
{
    EstimatorPyClass estimatorClass( m,
                                     "Estimator",
                                     R"doc(Class for orbit determination and covariance analysis.)doc" );

    bind_estimator_constructor_and_properties( estimatorClass );
    bind_estimator_execution_methods( estimatorClass );
}

}  // namespace estimation_analysis
}  // namespace estimation
}  // namespace tudatpy
