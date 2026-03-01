#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_estimation_analysis_ephemeris_fit_bridge.h"

#include <memory>
#include <vector>

#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/fitOrbitToEphemeris.h"

namespace py = pybind11;

namespace tba = tudat::basic_astrodynamics;
namespace tni = tudat::numerical_integrators;
namespace tss = tudat::simulation_setup;
namespace tep = tudat::estimatable_parameters;

namespace tudatpy
{
namespace estimation
{
namespace estimation_analysis
{

py::object create_best_fit_to_ephemeris_bridge(
    py::object bodies,
    py::object acceleration_models,
    py::object observed_bodies,
    py::object central_bodies,
    py::object integrator_settings,
    py::object initial_time,
    py::object final_time,
    py::object data_point_interval,
    py::object additional_parameter_names,
    py::object number_of_iterations,
    py::object reintegrate_variational_equations,
    py::object results_print_frequency )
{
    std::vector< std::shared_ptr< tep::EstimatableParameterSettings > > additionalParameters;
    if( !additional_parameter_names.is_none( ) )
    {
        additionalParameters = additional_parameter_names.cast<
            std::vector< std::shared_ptr< tep::EstimatableParameterSettings > > >( );
    }

    const auto accelerationModelMap = acceleration_models.cast< tba::AccelerationMap >( );
    const auto observedBodies = observed_bodies.cast< std::vector< std::string > >( );
    const auto centralBodies = central_bodies.cast< std::vector< std::string > >( );
    const auto integratorSettings = integrator_settings.cast< std::shared_ptr< tni::IntegratorSettings< TIME_TYPE > > >( );
    const auto initialTime = initial_time.cast< TIME_TYPE >( );
    const auto finalTime = final_time.cast< TIME_TYPE >( );
    const auto dataPointInterval = data_point_interval.cast< TIME_TYPE >( );
    const auto numberOfIterations = number_of_iterations.cast< int >( );
    const auto reintegrateVariationalEquations = reintegrate_variational_equations.cast< bool >( );
    const auto resultsPrintFrequency = results_print_frequency.cast< double >( );

    auto output = tss::createBestFitToCurrentEphemeris< TIME_TYPE, STATE_SCALAR_TYPE >(
        bodies.cast< const tss::SystemOfBodies& >( ),
        accelerationModelMap,
        observedBodies,
        centralBodies,
        integratorSettings,
        initialTime,
        finalTime,
        dataPointInterval,
        additionalParameters,
        numberOfIterations,
        reintegrateVariationalEquations,
        resultsPrintFrequency );

    return py::cast( output );
}

}  // namespace estimation_analysis
}  // namespace estimation
}  // namespace tudatpy
