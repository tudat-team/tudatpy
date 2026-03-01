#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_estimation_analysis_estimator_bindings.h"

#include <pybind11/pybind11.h>

#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "kernelExternTemplatesOrbitDetermination.h"

namespace py = pybind11;

namespace tudatpy
{
namespace estimation
{
namespace estimation_analysis
{

void bind_estimator_execution_methods( EstimatorPyClass& estimatorClass )
{
    estimatorClass
            .def( "perform_estimation",
                  &EstimatorType::estimateParameters,
                  py::arg( "estimation_input" ),
                  R"doc(Run a least-squares parameter estimation.)doc" )
            .def( "compute_covariance",
                  &EstimatorType::computeCovariance,
                  py::arg( "covariance_analysis_input" ),
                  R"doc(Compute the covariance matrix for the configured model and observations.)doc" );
}

}  // namespace estimation_analysis
}  // namespace estimation
}  // namespace tudatpy
