#ifndef TUDATPY_ESTIMATION_ANALYSIS_EPHEMERIS_FIT_BRIDGE_H
#define TUDATPY_ESTIMATION_ANALYSIS_EPHEMERIS_FIT_BRIDGE_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

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
    py::object results_print_frequency );

}  // namespace estimation_analysis
}  // namespace estimation
}  // namespace tudatpy

#endif
