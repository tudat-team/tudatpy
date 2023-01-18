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
#include "tudatpy/scalarTypes.h"

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



    py::class_<tp::PropagationPrintSettings,
            std::shared_ptr<tp::PropagationPrintSettings>>(m, "PropagationPrintSettings",
                                                           get_docstring("PropagationPrintSettings").c_str())
            .def_property("print_dependent_variable_indices",
                          &tp::PropagationPrintSettings::getPrintDependentVariableData,
                          &tp::PropagationPrintSettings::setPrintDependentVariableData,
                          get_docstring("PropagationPrintSettings.print_dependent_variable_indices").c_str() )
            .def_property("print_state_indices",
                          &tp::PropagationPrintSettings::getPrintStateData,
                          &tp::PropagationPrintSettings::setPrintStateData,
                          get_docstring("PropagationPrintSettings.print_state_indices").c_str() )
            .def_property("print_number_of_function_evaluations",
                          &tp::PropagationPrintSettings::getPrintNumberOfFunctionEvaluations,
                          &tp::PropagationPrintSettings::setPrintNumberOfFunctionEvaluations,
                          get_docstring("PropagationPrintSettings.print_number_of_function_evaluations").c_str() )
            .def_property("print_propagation_clock_time",
                          &tp::PropagationPrintSettings::getPrintPropagationTime,
                          &tp::PropagationPrintSettings::setPrintPropagationTime,
                          get_docstring("PropagationPrintSettings.print_propagation_clock_time").c_str() )
            .def_property("print_termination_reason",
                          &tp::PropagationPrintSettings::getPrintTerminationReason,
                          &tp::PropagationPrintSettings::setPrintTerminationReason,
                          get_docstring("PropagationPrintSettings.print_termination_reason").c_str() )
            .def_property("state_print_interval",
                          &tp::PropagationPrintSettings::getResultsPrintFrequencyInSteps,
                          &tp::PropagationPrintSettings::setResultsPrintFrequencyInSteps,
                          get_docstring("PropagationPrintSettings.state_print_interval").c_str() )
            .def_property("print_initial_and_final_conditions",
                          &tp::PropagationPrintSettings::getPrintInitialAndFinalConditions,
                          &tp::PropagationPrintSettings::setPrintInitialAndFinalConditions,
                          get_docstring("PropagationPrintSettings.print_initial_and_final_conditions").c_str() )
            .def("enable_all_printing",
                 py::overload_cast< >( &tp::PropagationPrintSettings::enableAllPrinting ),
                 get_docstring("PropagationPrintSettings.enable_all_printing").c_str() )
            .def("disable_all_printing",
                 &tp::PropagationPrintSettings::disableAllPrinting,
                 get_docstring("PropagationPrintSettings.disable_all_printing").c_str() );

    py::class_<tp::PropagatorProcessingSettings,
            std::shared_ptr<tp::PropagatorProcessingSettings>>(m, "PropagatorProcessingSettings",
                                                           get_docstring("PropagatorProcessingSettings").c_str())
            .def_property("set_integrated_result",
                          &tp::PropagatorProcessingSettings::getSetIntegratedResult,
                          &tp::PropagatorProcessingSettings::setIntegratedResult,
                          get_docstring("PropagatorProcessingSettings.set_integrated_result").c_str() )
            .def_property("clear_numerical_solution",
                          &tp::PropagatorProcessingSettings::getClearNumericalSolutions,
                          &tp::PropagatorProcessingSettings::setClearNumericalSolutions,
                          get_docstring("PropagatorProcessingSettings.clear_numerical_solution").c_str() );

    py::class_<tp::SingleArcPropagatorProcessingSettings,
            std::shared_ptr<tp::SingleArcPropagatorProcessingSettings>,
            tp::PropagatorProcessingSettings >(m, "SingleArcPropagatorProcessingSettings",
                                           get_docstring("SingleArcPropagatorProcessingSettings").c_str())
            .def_property_readonly("print_settings",
                                   &tp::SingleArcPropagatorProcessingSettings::getPrintSettings,
                                   get_docstring("SingleArcPropagatorProcessingSettings.print_settings").c_str() )
            .def_property("results_save_frequency_in_steps",
                          &tp::SingleArcPropagatorProcessingSettings::getResultsSaveFrequencyInSteps,
                          &tp::SingleArcPropagatorProcessingSettings::setResultsSaveFrequencyInSteps,
                          get_docstring("SingleArcPropagatorProcessingSettings.results_save_frequency_in_steps").c_str() )
            .def_property("results_save_frequency_in_seconds",
                          &tp::SingleArcPropagatorProcessingSettings::getResultsSaveFrequencyInSeconds,
                          &tp::SingleArcPropagatorProcessingSettings::setResultsSaveFrequencyInSeconds,
                          get_docstring("SingleArcPropagatorProcessingSettings.results_save_frequency_in_seconds").c_str() );

    py::class_<tp::MultiArcPropagatorProcessingSettings,
            std::shared_ptr<tp::MultiArcPropagatorProcessingSettings>,
            tp::PropagatorProcessingSettings >(m, "MultiArcPropagatorProcessingSettings",
                                           get_docstring("MultiArcPropagatorProcessingSettings").c_str())
            .def("set_print_settings_for_all_arcs",
                 &tp::MultiArcPropagatorProcessingSettings::resetAndApplyConsistentSingleArcPrintSettings,
                 get_docstring("MultiArcPropagatorProcessingSettings.disable_all_printing").c_str() )
            .def_property("print_output_on_first_arc_only",
                 &tp::MultiArcPropagatorProcessingSettings::getPrintFirstArcOnly,
                 &tp::MultiArcPropagatorProcessingSettings::resetPrintFirstArcOnly,
                 get_docstring("MultiArcPropagatorProcessingSettings.print_output_on_first_arc_only").c_str() )
            .def_property_readonly("single_arc_settings",
                          &tp::MultiArcPropagatorProcessingSettings::getSingleArcSettings,
                          get_docstring("MultiArcPropagatorProcessingSettings.single_arc_settings").c_str() );


    py::class_<tp::HybridArcPropagatorProcessingSettings,
            std::shared_ptr<tp::HybridArcPropagatorProcessingSettings>,
            tp::PropagatorProcessingSettings >(m, "HybridArcPropagatorProcessingSettings",
                                           get_docstring("HybridArcPropagatorProcessingSettings").c_str());

    ///////////////////////////////////////////////////////////////////////////////////////

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
    py::class_<
            tp::PropagatorSettings<double>,
            std::shared_ptr<tp::PropagatorSettings<double>>>(m, "PropagatorSettings",
                                                             get_docstring("PropagatorSettings").c_str())
            .def_property("initial_states",
                          &tp::PropagatorSettings<double>::getInitialStates,
                          &tp::PropagatorSettings<double>::resetInitialStates,
                          get_docstring("PropagatorSettings.initial_states").c_str());

    py::class_<
            tp::MultiArcPropagatorSettings<double,TIME_TYPE>,
            std::shared_ptr<tp::MultiArcPropagatorSettings<double,TIME_TYPE>>,
            tp::PropagatorSettings<double>>(m, "MultiArcPropagatorSettings",
                                            get_docstring("MultiArcPropagatorSettings").c_str());

    py::class_<
            tp::HybridArcPropagatorSettings<double,TIME_TYPE>,
            std::shared_ptr<tp::HybridArcPropagatorSettings<double,TIME_TYPE>>,
            tp::PropagatorSettings<double>>(m, "HybridArcPropagatorSettings",
                                            get_docstring("HybridArcPropagatorSettings").c_str());

    py::class_<
            tp::SingleArcPropagatorSettings<double,TIME_TYPE>,
            std::shared_ptr<tp::SingleArcPropagatorSettings<double,TIME_TYPE>>,
            tp::PropagatorSettings<double>>(m, "SingleArcPropagatorSettings",
                                            get_docstring("SingleArcPropagatorSettings").c_str())
            .def_property("termination_settings",
                          &tp::SingleArcPropagatorSettings<double,TIME_TYPE>::getTerminationSettings,
                          &tp::SingleArcPropagatorSettings<double,TIME_TYPE>::resetTerminationSettings,
                          get_docstring("SingleArcPropagatorSettings.termination_settings").c_str() )
            .def_property_readonly("processing_settings",
                                   &tp::SingleArcPropagatorSettings<double,TIME_TYPE>::getOutputSettings )
            .def_property_readonly("print_settings",
                                   &tp::SingleArcPropagatorSettings<double,TIME_TYPE>::getPrintSettings );


    py::class_<
            tp::TranslationalStatePropagatorSettings<double,TIME_TYPE>,
            std::shared_ptr<tp::TranslationalStatePropagatorSettings<double,TIME_TYPE>>,
            tp::SingleArcPropagatorSettings<double, TIME_TYPE>>(m, "TranslationalStatePropagatorSettings",
                                                     get_docstring("TranslationalStatePropagatorSettings").c_str())

            .def("get_propagated_state_size",
                 &tp::TranslationalStatePropagatorSettings<double,TIME_TYPE>::getPropagatedStateSize)
            .def("reset_and_recreate_acceleration_models",
                 &tp::TranslationalStatePropagatorSettings<double,TIME_TYPE>::resetAccelerationModelsMap,
                 py::arg("new_acceleration_settings"),
                 py::arg("bodies"));


    //            .def_property_readonly("acceleration_settings",
    //                                   &tp::TranslationalStatePropagatorSettings<double>::getAccelerationSettingsMap);


    py::class_<
            tp::MultiTypePropagatorSettings<double,TIME_TYPE>,
            std::shared_ptr<tp::MultiTypePropagatorSettings<double,TIME_TYPE>>,
            tp::SingleArcPropagatorSettings<double,TIME_TYPE>>(m, "MultiTypePropagatorSettings",
                                                     get_docstring("MultiTypePropagatorSettings").c_str())
            .def("reset_initial_states", &tp::MultiTypePropagatorSettings<double,TIME_TYPE>::resetInitialStates,
                 py::arg("initial_states"))
            .def("recreate_state_derivative_models",
                 &tp::MultiTypePropagatorSettings<double,TIME_TYPE>::resetIntegratedStateModels,
                 py::arg("bodies"))
            .def("single_type_settings", &tp::MultiTypePropagatorSettings<double,TIME_TYPE>::getSingleTypePropagatorSettings,
                 py::arg("state_type"))
            .def_property_readonly("propagator_settings_per_type",
                                   &tp::MultiTypePropagatorSettings<double,TIME_TYPE>::getPropagatorSettingsMap,
                                   get_docstring("MultiTypePropagatorSettings.propagator_settings_per_type").c_str());

    py::class_<
            tp::RotationalStatePropagatorSettings<double,TIME_TYPE>,
            std::shared_ptr<tp::RotationalStatePropagatorSettings<double,TIME_TYPE>>,
            tp::SingleArcPropagatorSettings<double,TIME_TYPE>>(m, "RotationalStatePropagatorSettings",
                                                     get_docstring("RotationalStatePropagatorSettings").c_str());

    py::class_<
            tp::MassPropagatorSettings<double,TIME_TYPE>,
            std::shared_ptr<tp::MassPropagatorSettings<double,TIME_TYPE>>,
            tp::SingleArcPropagatorSettings<double,TIME_TYPE>>(m, "MassPropagatorSettings",
                                                     get_docstring("MassPropagatorSettings").c_str());

    py::class_<
            tp::CustomStatePropagatorSettings<double, TIME_TYPE>,
            std::shared_ptr<tp::CustomStatePropagatorSettings<double, TIME_TYPE>>,
            tp::SingleArcPropagatorSettings<double,TIME_TYPE>>(m, "CustomStatePropagatorSettings",
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
          const double>(&tp::translationalStatePropagatorSettingsDeprecated<double,TIME_TYPE>),
          py::arg("central_bodies"),
          py::arg("acceleration_models"),
          py::arg("bodies_to_integrate"),
          py::arg("initial_states"),
          py::arg("termination_settings"),
          py::arg("propagator") = tp::cowell,
          py::arg("output_variables") = std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >(),
          py::arg("print_interval") = TUDAT_NAN );

    m.def("translational",
          tp::translationalStatePropagatorSettings<double,TIME_TYPE>,
          py::arg("central_bodies"),
          py::arg("acceleration_models"),
          py::arg("bodies_to_integrate"),
          py::arg("initial_states"),
          py::arg("initial_time"),
          py::arg("integrator_settings"),
          py::arg("termination_settings"),
          py::arg("propagator") = tp::cowell,
          py::arg("output_variables") = std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >(),
          py::arg("processing_settings") = std::make_shared< tp::SingleArcPropagatorProcessingSettings >( ),
          get_docstring("translational").c_str());

    m.def("rotational",
          py::overload_cast<
          const tba::TorqueModelMap &,
          const std::vector<std::string> &,
          const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
          const std::shared_ptr<tp::PropagationTerminationSettings>,
          const tp::RotationalPropagatorType,
          const std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >,
          const double>(&tp::rotationalStatePropagatorSettingsDeprecated<double,TIME_TYPE>),
          py::arg("torque_models"),
          py::arg("bodies_to_integrate"),
          py::arg("initial_states"),
          py::arg("termination_settings"),
          py::arg("propagator") = tp::quaternions,
          py::arg("output_variables") = std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >(),
          py::arg("print_interval") = TUDAT_NAN );

    m.def("rotational",
          &tp::rotationalStatePropagatorSettings<double,TIME_TYPE>,
          py::arg("torque_models"),
          py::arg("bodies_to_integrate"),
          py::arg("initial_states"),
          py::arg("initial_time"),
          py::arg("integrator_settings"),
          py::arg("termination_settings"),
          py::arg("propagator") = tp::quaternions,
          py::arg("output_variables") = std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >(),
          py::arg("processing_settings") = std::make_shared< tp::SingleArcPropagatorProcessingSettings >( ),
          get_docstring("rotational").c_str());

    m.def("mass",
          py::overload_cast<
          const std::vector<std::string>,
          const std::map<std::string, std::vector<std::shared_ptr<tba::MassRateModel> > > &,
          const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
          const std::shared_ptr<tp::PropagationTerminationSettings>,
          const std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >,
          const double>(&tp::massPropagatorSettingsDeprecated<double,TIME_TYPE>),
          py::arg("bodies_with_mass_to_propagate"),
          py::arg("mass_rate_models"),
          py::arg("initial_body_masses"),
          py::arg("termination_settings"),
          py::arg("output_variables") = std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >(),
          py::arg("print_interval") = TUDAT_NAN );

    m.def("mass",
          &tp::massPropagatorSettings<double,TIME_TYPE>,
          py::arg("bodies_with_mass_to_propagate"),
          py::arg("mass_rate_models"),
          py::arg("initial_body_masses"),
          py::arg("initial_time"),
          py::arg("integrator_settings"),
          py::arg("termination_settings"),
          py::arg("output_variables") = std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >(),
          py::arg("processing_settings") = std::make_shared< tp::SingleArcPropagatorProcessingSettings >( ),
          get_docstring("mass").c_str());

    m.def("custom_state",
          &tp::customStatePropagatorSettings<double, TIME_TYPE>,
          py::arg("state_derivative_function"),
          py::arg("initial_state"),
          py::arg("initial_time"),
          py::arg("integrator_settings"),
          py::arg("termination_settings"),
          py::arg("output_variables") = std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >(),
          py::arg("processing_settings") = std::make_shared< tp::SingleArcPropagatorProcessingSettings >( ),
          get_docstring("custom_state").c_str());




    m.def("multitype",
          py::overload_cast<
          const std::vector<std::shared_ptr<tp::SingleArcPropagatorSettings<double,TIME_TYPE> > >,
          const std::shared_ptr<tp::PropagationTerminationSettings>,
          const std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >,
          const double>(&tp::multiTypePropagatorSettingsDeprecated<double,TIME_TYPE>),
          py::arg("propagator_settings_list"),
          py::arg("termination_settings"),
          py::arg("output_variables") = std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >(),
          py::arg("print_interval") = TUDAT_NAN );

    m.def("multitype",
          &tp::multiTypePropagatorSettings<double,TIME_TYPE>,
          py::arg("propagator_settings_list"),
          py::arg("integrator_settings"),
          py::arg("initial_time"),
          py::arg("termination_settings"),
          py::arg("output_variables") = std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >(),
          py::arg("processing_settings") = std::make_shared< tp::SingleArcPropagatorProcessingSettings >( ),
          get_docstring("multitype").c_str());


    m.def("multi_arc",
          &tp::multiArcPropagatorSettings<double,TIME_TYPE>,
          py::arg("single_arc_settings"),
          py::arg("transfer_state_to_next_arc") = false,
          py::arg("processing_settings") = std::make_shared< tp::MultiArcPropagatorProcessingSettings >( ),
          get_docstring("multi_arc").c_str());

    m.def("hybrid_arc",
          &tp::hybridArcPropagatorSettings<double,TIME_TYPE>,
          py::arg("single_arc_settings"),
          py::arg("multi_arc_settings"),
          py::arg("processing_settings") = std::make_shared< tp::HybridArcPropagatorProcessingSettings >( ),
          get_docstring("hybrid_arc").c_str());

    ///////////////////////////////////////////////////////////////////////////////////////


    py::class_<tp::PropagationTerminationSettings,
            std::shared_ptr<tp::PropagationTerminationSettings>>(m, "PropagationTerminationSettings",
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


