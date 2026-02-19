#ifndef TUDATPY_EXPOSE_ESTIMATION_ANALYSIS_ESTIMATOR_BINDINGS_H
#define TUDATPY_EXPOSE_ESTIMATION_ANALYSIS_ESTIMATOR_BINDINGS_H

#include <memory>

#include <pybind11/pybind11.h>

#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/estimationInterfacesForwardDeclarations.h"

namespace py = pybind11;

namespace tudatpy
{
namespace estimation
{
namespace estimation_analysis
{

using EstimatorType = tudat::simulation_setup::OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE, 0 >;
using EstimatorPyClass = py::class_< EstimatorType, std::shared_ptr< EstimatorType > >;

void bind_estimator_constructor_and_properties( EstimatorPyClass& estimatorClass );
void bind_estimator_execution_methods( EstimatorPyClass& estimatorClass );

}  // namespace estimation_analysis
}  // namespace estimation
}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_ESTIMATION_ANALYSIS_ESTIMATOR_BINDINGS_H
