/*    Copyright (c) 2010-2021, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_propagator_setup.h"

#include "tudatpy/docstrings.h"

#include <tudat/simulation/propagation_setup.h>
#include <tudat/astro/propagators/getZeroProperModeRotationalInitialState.h>

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
namespace tba = tudat::basic_astrodynamics;
namespace tss = tudat::simulation_setup;
namespace tp = tudat::propagators;
namespace tinterp = tudat::interpolators;
namespace te = tudat::ephemerides;
namespace tni = tudat::numerical_integrators;
namespace trf = tudat::reference_frames;
namespace tmrf = tudat::root_finders;


namespace tudatpy {
namespace numerical_simulation {
namespace propagation_setup {
namespace propagator {

void expose_propagator_setup(py::module &m) {

    // ENUMS
    py::enum_<tp::TranslationalPropagatorType>(m, "TranslationalPropagatorType",
                                               get_docstring("TranslationalPropagatorType").c_str())
            .value("undefined_translational_propagator",
                   tp::TranslationalPropagatorType::undefined_translational_propagator,
                   get_docstring("TranslationalPropagatorType.undefined_translational_propagator").c_str())
            .value("cowell",
                   tp::TranslationalPropagatorType::cowell,
                   get_docstring("TranslationalPropagatorType.cowell").c_str())
            .value("encke",
                   tp::TranslationalPropagatorType::encke,
                   get_docstring("TranslationalPropagatorType.encke").c_str())
            .value("gauss_keplerian",
                   tp::TranslationalPropagatorType::gauss_keplerian,
                   get_docstring("TranslationalPropagatorType.gauss_keplerian").c_str())
            .value("gauss_modified_equinoctial",
                   tp::TranslationalPropagatorType::gauss_modified_equinoctial,
                   get_docstring("TranslationalPropagatorType.gauss_modified_equinoctial").c_str())
            .value("unified_state_model_quaternions",
                   tp::TranslationalPropagatorType::unified_state_model_quaternions,
                   get_docstring("TranslationalPropagatorType.unified_state_model_quaternions").c_str())
            .value("unified_state_model_modified_rodrigues_parameters",
                   tp::TranslationalPropagatorType::unified_state_model_modified_rodrigues_parameters,
                   get_docstring("TranslationalPropagatorType.unified_state_model_modified_rodrigues_parameters").c_str())
            .value("unified_state_model_exponential_map",
                   tp::unified_state_model_exponential_map,
                   get_docstring("TranslationalPropagatorType.unified_state_model_exponential_map").c_str())
            .export_values();

    py::enum_<tp::RotationalPropagatorType>(m, "RotationalPropagatorType",
                                            get_docstring("RotationalPropagatorType").c_str())
            .value("undefined_rotational_propagator",
                   tp::RotationalPropagatorType::undefined_rotational_propagator,
                   get_docstring("RotationalPropagatorType.undefined_rotational_propagator").c_str())
            .value("quaternions",
                   tp::RotationalPropagatorType::quaternions,
                   get_docstring("RotationalPropagatorType.quaternions").c_str())
            .value("modified_rodrigues_parameters",
                   tp::RotationalPropagatorType::modified_rodrigues_parameters,
                   get_docstring("RotationalPropagatorType.modified_rodrigues_parameters").c_str())
            .value("exponential_map",
                   tp::RotationalPropagatorType::exponential_map,
                   get_docstring("RotationalPropagatorType.exponential_map").c_str())
            .export_values();

    py::enum_<tp::PropagationTerminationTypes>(m, "PropagationTerminationTypes",
                                               get_docstring("PropagationTerminationTypes").c_str())
            .value("time_stopping_condition_type",
                   tp::PropagationTerminationTypes::time_stopping_condition,
                   get_docstring("PropagationTerminationTypes.time_stopping_condition_type").c_str())
            .value("cpu_time_stopping_condition_type",
                   tp::PropagationTerminationTypes::cpu_time_stopping_condition,
                   get_docstring("PropagationTerminationTypes.cpu_time_stopping_condition_type").c_str())
            .value("dependent_variable_stopping_condition_type",
                   tp::PropagationTerminationTypes::dependent_variable_stopping_condition,
                   get_docstring("PropagationTerminationTypes.dependent_variable_stopping_condition_type").c_str())
            .value("hybrid_stopping_condition_type",
                   tp::PropagationTerminationTypes::hybrid_stopping_condition,
                   get_docstring("PropagationTerminationTypes.hybrid_stopping_condition_type").c_str())
            .value("custom_stopping_condition_type",
                   tp::PropagationTerminationTypes::custom_stopping_condition,
                   get_docstring("PropagationTerminationTypes.custom_stopping_condition_type").c_str())
            .export_values();

    py::enum_<tp::IntegratedStateType>(m, "StateType",
                                       get_docstring("StateType").c_str())
            .value("hybrid_type", tp::IntegratedStateType::hybrid,
                   get_docstring("StateType.hybrid_type").c_str())
            .value("translational_type", tp::IntegratedStateType::translational_state,
                   get_docstring("StateType.translational_type").c_str())
            .value("rotational_type", tp::IntegratedStateType::rotational_state,
                   get_docstring("StateType.rotational_type").c_str())
            .value("mass_type", tp::IntegratedStateType::body_mass_state,
                   get_docstring("StateType.mass_type").c_str())
            .value("custom_type", tp::IntegratedStateType::custom_state,
                   get_docstring("StateType.custom_type").c_str())
            .export_values();


    //        py::enum_<tp::VariableType>(m, "VariableType")
    //                .value("independent_variable", tp::VariableType::independentVariable)
    //                .value("cpu_time_variable", tp::VariableType::cpuTimeVariable)
    //                .value("state_variable", tp::VariableType::stateVariable)
    //                .value("dependent_variable", tp::VariableType::dependentVariable)
    //                .export_values();


    // CLASSES
    py::class_<tp::DependentVariableSaveSettings,
            std::shared_ptr<tp::DependentVariableSaveSettings>>(m, "DependentVariableSaveSettings",
                                                                get_docstring("DependentVariableSaveSettings").c_str());
    //                .def(py::init<
    //                             const std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings>>,
    //                             const bool>(),
    //                     py::arg("dependent_variables"),
    //                     py::arg("print_dependent_variable_types") = true);

    py::class_<
            tp::PropagatorSettings<double>,
            std::shared_ptr<tp::PropagatorSettings<double>>>(m, "PropagatorSettings",
                                                             get_docstring("PropagatorSettings").c_str())
            .def_property("initial_states",
                          &tp::PropagatorSettings<double>::getInitialStates,
                          &tp::PropagatorSettings<double>::resetInitialStates,
                          get_docstring("PropagatorSettings.initial_states").c_str());

    py::class_<
            tp::MultiArcPropagatorSettings<double>,
            std::shared_ptr<tp::MultiArcPropagatorSettings<double>>,
            tp::PropagatorSettings<double>>(m, "MultiArcPropagatorSettings",
                                            get_docstring("MultiArcPropagatorSettings").c_str());

    py::class_<
            tp::HybridArcPropagatorSettings<double>,
            std::shared_ptr<tp::HybridArcPropagatorSettings<double>>,
            tp::PropagatorSettings<double>>(m, "HybridArcPropagatorSettings",
                    get_docstring("HybridArcPropagatorSettings").c_str());

    py::class_<
            tp::SingleArcPropagatorSettings<double>,
            std::shared_ptr<tp::SingleArcPropagatorSettings<double>>,
            tp::PropagatorSettings<double>>(m, "SingleArcPropagatorSettings",
                                            get_docstring("SingleArcPropagatorSettings").c_str())
            .def_property("termination_settings",
                          &tp::SingleArcPropagatorSettings<double>::getTerminationSettings,
                          &tp::SingleArcPropagatorSettings<double>::resetTerminationSettings,
                          get_docstring("SingleArcPropagatorSettings.termination_settings").c_str());

    py::class_<
            tp::TranslationalStatePropagatorSettings<double>,
            std::shared_ptr<tp::TranslationalStatePropagatorSettings<double>>,
            tp::SingleArcPropagatorSettings<double>>(m, "TranslationalStatePropagatorSettings",
                                                     get_docstring("TranslationalStatePropagatorSettings").c_str())
            //                .def(// ctor 1
            //                        py::init<
            //                                const std::vector<std::string> &,
            //                                const tba::AccelerationMap &,
            //                                const std::vector<std::string> &,
            //                                const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
            //                                const std::shared_ptr<tp::PropagationTerminationSettings>,
            //                                const tp::TranslationalPropagatorType,
            //                                const std::shared_ptr<tp::DependentVariableSaveSettings>,
            //                                const double>(),
            //                        py::arg("central_bodies"),
            //                        py::arg("acceleration_models"),
            //                        py::arg("bodies_to_integrate"),
            //                        py::arg("initial_states"),
            //                        py::arg("termination_settings"),
            //                        py::arg("propagator") = tp::TranslationalPropagatorType::cowell,
            //                        py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
            //                        py::arg("print_interval") = TUDAT_NAN)
            //                .def(// ctor 2
            //                        py::init<const std::vector<std::string> &,
            //                                const tss::SelectedAccelerationMap &,
            //                                const std::vector<std::string> &,
            //                                const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
            //                                const std::shared_ptr<tp::PropagationTerminationSettings>,
            //                                const tp::TranslationalPropagatorType,
            //                                const std::shared_ptr<tp::DependentVariableSaveSettings>,
            //                                const double>(),
            //                        py::arg("central_bodies"),
            //                        py::arg("acceleration_settings"),
            //                        py::arg("bodies_to_integrate"),
            //                        py::arg("initial_states"),
            //                        py::arg("termination_settings"),
            //                        py::arg("propagator") = tp::cowell,
            //                        py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
            //                        py::arg("print_interval") = TUDAT_NAN)
            //                .def(// ctor 3
            //                        py::init<const std::vector<std::string> &,
            //                                const tba::AccelerationMap &,
            //                                const std::vector<std::string> &,
            //                                const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
            //                                const double,
            //                                const tp::TranslationalPropagatorType,
            //                                const std::shared_ptr<tp::DependentVariableSaveSettings>,
            //                                const double>(),
            //                        py::arg("central_bodies"),
            //                        py::arg("acceleration_models"),
            //                        py::arg("bodies_to_integrate"),
            //                        py::arg("initial_states"),
            //                        py::arg("termination_time"),
            //                        py::arg("propagator") = tp::cowell,
            //                        py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
            //                        py::arg("print_interval") = TUDAT_NAN)
            //                .def(// ctor 4
            //                        py::init<const std::vector<std::string> &,
            //                                const tss::SelectedAccelerationMap &,
            //                                const std::vector<std::string> &,
            //                                const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
            //                                const double,
            //                                const tp::TranslationalPropagatorType,
            //                                const std::shared_ptr<tp::DependentVariableSaveSettings>,
            //                                const double>(),
            //                        py::arg("central_bodies"),
            //                        py::arg("acceleration_settings"),
            //                        py::arg("bodies_to_integrate"),
            //                        py::arg("initial_states"),
            //                        py::arg("termination_time"),
            //                        py::arg("propagator") = tp::cowell,
            //                        py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
            //                        py::arg("print_interval") = TUDAT_NAN)
//            .def("recreate_state_derivative_models",
//                 &tp::TranslationalStatePropagatorSettings<double>::resetIntegratedStateModels,
//                 py::arg("bodies"))
            .def("get_propagated_state_size",
                 &tp::TranslationalStatePropagatorSettings<double>::getPropagatedStateSize)
            .def("reset_and_recreate_acceleration_models",
                 &tp::TranslationalStatePropagatorSettings<double>::resetAccelerationModelsMap,
                 py::arg("new_acceleration_settings"),
                 py::arg("bodies"));
//            .def_property_readonly("acceleration_settings",
//                                   &tp::TranslationalStatePropagatorSettings<double>::getAccelerationSettingsMap);


    py::class_<
            tp::MultiTypePropagatorSettings<double>,
            std::shared_ptr<tp::MultiTypePropagatorSettings<double>>,
            tp::SingleArcPropagatorSettings<double>>(m, "MultiTypePropagatorSettings",
                                                     get_docstring("MultiTypePropagatorSettings").c_str())
            .def("reset_initial_states", &tp::MultiTypePropagatorSettings<double>::resetInitialStates,
                 py::arg("initial_states"))
            .def("recreate_state_derivative_models",
                 &tp::MultiTypePropagatorSettings<double>::resetIntegratedStateModels,
                 py::arg("bodies"))
            .def("single_type_settings", &tp::MultiTypePropagatorSettings<double>::getSingleTypePropagatorSettings,
                 py::arg("state_type"))
            .def_property_readonly("propagator_settings_per_type",
                                   &tp::MultiTypePropagatorSettings<double>::getPropagatorSettingsMap,
                                   get_docstring("MultiTypePropagatorSettings.propagator_settings_per_type").c_str());

    py::class_<
            tp::RotationalStatePropagatorSettings<double>,
            std::shared_ptr<tp::RotationalStatePropagatorSettings<double>>,
            tp::SingleArcPropagatorSettings<double>>(m, "RotationalStatePropagatorSettings",
                                                     get_docstring("RotationalStatePropagatorSettings").c_str());

    py::class_<
            tp::MassPropagatorSettings<double>,
            std::shared_ptr<tp::MassPropagatorSettings<double>>,
            tp::SingleArcPropagatorSettings<double>>(m, "MassPropagatorSettings",
                                                     get_docstring("MassPropagatorSettings").c_str());

    py::class_<
            tp::CustomStatePropagatorSettings<double, double>,
            std::shared_ptr<tp::CustomStatePropagatorSettings<double, double>>,
            tp::SingleArcPropagatorSettings<double>>(m, "CustomStatePropagatorSettings",
                                                     get_docstring("CustomStatePropagatorSettings").c_str());



    m.def("translational",
          py::overload_cast<
          const std::vector<std::string> &,
          const tba::AccelerationMap &,
          const std::vector<std::string> &,
          const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
          const std::shared_ptr<tp::PropagationTerminationSettings>,
          const tp::TranslationalPropagatorType,
          const std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> > &,
          const double>(&tp::translationalStatePropagatorSettings<double>),
          py::arg("central_bodies"),
          py::arg("acceleration_models"),
          py::arg("bodies_to_integrate"),
          py::arg("initial_states"),
          py::arg("termination_settings"),
          py::arg("propagator") = tp::cowell,
          py::arg("output_variables") = std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >(),
          py::arg("print_interval") = TUDAT_NAN,
          get_docstring("translational").c_str());

    m.def("rotational",
          py::overload_cast<
          const tba::TorqueModelMap &,
          const std::vector<std::string> &,
          const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
          const std::shared_ptr<tp::PropagationTerminationSettings>,
          const tp::RotationalPropagatorType,
          const std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >,
          const double>(&tp::rotationalStatePropagatorSettings<double>),
          py::arg("torque_models"),
          py::arg("bodies_to_integrate"),
          py::arg("initial_states"),
          py::arg("termination_settings"),
          py::arg("propagator") = tp::quaternions,
          py::arg("output_variables") = std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >(),
          py::arg("print_interval") = TUDAT_NAN,
          get_docstring("rotational").c_str());


    m.def("mass",
          py::overload_cast<
          const std::vector<std::string>,
          const std::map<std::string, std::vector<std::shared_ptr<tba::MassRateModel> > > &,
          const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
          const std::shared_ptr<tp::PropagationTerminationSettings>,
          const std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >,
          const double>(&tp::massPropagatorSettings<double>),
          py::arg("bodies_with_mass_to_propagate"),
          py::arg("mass_rate_models"),
          py::arg("initial_body_masses"),
          py::arg("termination_settings"),
          py::arg("output_variables") = std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >(),
          py::arg("print_interval") = TUDAT_NAN,
          get_docstring("mass").c_str());

    m.def("custom_state",
          &tp::customStatePropagatorSettings<double, double>,
          py::arg("state_derivative_function"),
          py::arg("initial_state"),
          py::arg("termination_settings"),
          py::arg("output_variables") = std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >(),
          py::arg("print_interval") = TUDAT_NAN,
          get_docstring("custom_state").c_str());




    m.def("multitype",
          py::overload_cast<
          const std::vector<std::shared_ptr<tp::SingleArcPropagatorSettings<double> > >,
          const std::shared_ptr<tp::PropagationTerminationSettings>,
          const std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >,
          const double>(&tp::multiTypePropagatorSettings<double>),
          py::arg("propagator_settings_list"),
          py::arg("termination_settings"),
          py::arg("output_variables") = std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >(),
          py::arg("print_interval") = TUDAT_NAN,
          get_docstring("multitype").c_str());


    m.def("multi_arc",
          &tp::multiArcPropagatorSettings<double>,
          py::arg("single_arc_settings"),
          py::arg("transfer_state_to_next_arc") = false,
          get_docstring("multi_arc").c_str());

    m.def("hybrid_arc",
          &tp::hybridArcPropagatorSettings<double>,
          py::arg("single_arc_settings"),
          py::arg("multi_arc_settings"),
          get_docstring("hybrid_arc").c_str());

    py::class_<tp::PropagationTerminationSettings,
            std::shared_ptr<tp::PropagationTerminationSettings>>
            PropagationTerminationSettings_(m, "PropagationTerminationSettings",
                                            get_docstring("PropagationTerminationSettings").c_str());

    py::class_<
            tp::PropagationDependentVariableTerminationSettings,
            std::shared_ptr<tp::PropagationDependentVariableTerminationSettings>,
            tp::PropagationTerminationSettings>(m, "PropagationDependentVariableTerminationSettings",
                                                get_docstring(
                                                    "PropagationDependentVariableTerminationSettings").c_str());

    py::class_<
            tp::PropagationTimeTerminationSettings,
            std::shared_ptr<tp::PropagationTimeTerminationSettings>,
            tp::PropagationTerminationSettings>(m, "PropagationTimeTerminationSettings",
                                                get_docstring(
                                                    "PropagationTimeTerminationSettings").c_str());

    py::class_<
            tp::PropagationCPUTimeTerminationSettings,
            std::shared_ptr<tp::PropagationCPUTimeTerminationSettings>,
            tp::PropagationTerminationSettings>(m, "PropagationCPUTimeTerminationSettings",
                                                get_docstring(
                                                    "PropagationCPUTimeTerminationSettings").c_str());

    py::class_<
            tp::PropagationCustomTerminationSettings,
            std::shared_ptr<tp::PropagationCustomTerminationSettings>,
            tp::PropagationTerminationSettings>(m, "PropagationCustomTerminationSettings",
                                                get_docstring(
                                                    "PropagationCustomTerminationSettings").c_str());

    py::class_<
            tp::PropagationHybridTerminationSettings,
            std::shared_ptr<tp::PropagationHybridTerminationSettings>,
            tp::PropagationTerminationSettings>(m, "PropagationHybridTerminationSettings",
                                                get_docstring(
                                                    "PropagationHybridTerminationSettings").c_str());

    //                .def(py::init<
    //                             const std::shared_ptr<tp::SingleDependentVariableSaveSettings>,
    //                             const double,
    //                             const bool,
    //                             const bool,
    //                             const std::shared_ptr<tudat::root_finders::RootFinderSettings>>(),
    //                     py::arg("dependent_variadble_settings"),
    //                     py::arg("limit_value"),
    //                     py::arg("use_as_lower_limit"),
    //                     py::arg("terminate_exactly_on_final_condition") = false,
    //                     py::arg("termination_root_finder_settings") = nullptr);

    m.def("time_termination",
          &tp::propagationTimeTerminationSettings,
          py::arg("termination_time"),
          py::arg("terminate_exactly_on_final_condition") = false,
          get_docstring("time_termination").c_str());

    m.def("cpu_time_termination",
          &tp::propagationCPUTimeTerminationSettings,
          py::arg("cpu_termination_time"),
          get_docstring("cpu_time_termination").c_str());

    m.def("dependent_variable_termination",
          &tp::propagationDependentVariableTerminationSettings,
          py::arg("dependent_variable_settings"),
          py::arg("limit_value"),
          py::arg("use_as_lower_limit"),
          py::arg("terminate_exactly_on_final_condition") = false,
          py::arg("termination_root_finder_settings") = nullptr,
          get_docstring("dependent_variable_termination").c_str());

    m.def("custom_termination",
          &tp::popagationCustomTerminationSettings,
          py::arg("custom_condition"),
          get_docstring("custom_termination").c_str());

    m.def("hybrid_termination",
          &tp::propagationHybridTerminationSettings,
          py::arg("termination_settings"),
          py::arg("fulfill_single_condition"),
          get_docstring("hybrid_termination").c_str());

    m.def("add_dependent_variable_settings",
          &tp::addDepedentVariableSettings< double >,
          py::arg("dependent_variable_settings"),
          py::arg("propagator_settings"),
          get_docstring("add_dependent_variable_settings").c_str());


}

}// namespace propagator
}// namespace propagation_setup
}// namespace numerical_simulation
}// namespace tudatpy


