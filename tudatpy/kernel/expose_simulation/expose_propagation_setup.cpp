/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_propagation_setup.h"

#include "expose_propagation_setup/expose_dependent_variable_setup.h"
#include "expose_propagation_setup/expose_torque_setup.h"
#include "expose_propagation_setup/expose_acceleration_setup.h"

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


void expose_mass_rate_setup(py::module &m)
{
    py::enum_<tba::AvailableMassRateModels>(m, "AvailableMassRateModels")
            .value("undefined_mass_rate_type", tba::AvailableMassRateModels::undefined_mass_rate_model)
            .value("custom_mass_rate_type", tba::AvailableMassRateModels::custom_mass_rate_model)
            .value("from_thrust_mass_rate_type", tba::AvailableMassRateModels::from_thrust_mass_rate_model)
            .export_values();

    py::class_<tss::MassRateModelSettings,
            std::shared_ptr<tss::MassRateModelSettings>>(m, "MassRateModelSettings")
            .def(py::init<const tudat::basic_astrodynamics::AvailableMassRateModels>(),
                 py::arg("mass_rate_type"));

    py::class_<tss::FromThrustMassModelSettings,
            std::shared_ptr<tss::FromThrustMassModelSettings>,
            tss::MassRateModelSettings>(m, "FromThrustMassModelSettings")
            .def(py::init<const bool, const std::string&>(),
                 py::arg("use_all_thrust_models") = 1,
                 py::arg("associated_thrust_source") = "" );

    m.def("custom", &tss::customMassRate,
          py::arg( "mass_rate_function" ) );

    m.def("from_thrust", &tss::fromThrustMassRate,
          py::arg( "use_all_thrust_models" ) = 1,
          py::arg( "associated_thrust_source" ) = "" );

}

void expose_integrator_setup(py::module &m) {


    py::enum_<tni::AvailableIntegrators>(m, "AvailableIntegrators")
            .value("euler_type", tni::AvailableIntegrators::euler)
            .value("runge_kutta_4_type", tni::AvailableIntegrators::rungeKutta4)
            .value("runge_kutta_variable_step_size_type", tni::AvailableIntegrators::rungeKuttaVariableStepSize)
            .value("bulirsch_stoer_type", tni::AvailableIntegrators::bulirschStoer)
            .value("adams_bashforth_moulton_type", tni::AvailableIntegrators::adamsBashforthMoulton)
            .export_values();

    py::enum_<tni::RungeKuttaCoefficients::CoefficientSets>(m, "RKCoefficientSets")
            .value("rkf_45", tni::RungeKuttaCoefficients::rungeKuttaFehlberg45)
            .value("rkf_56", tni::RungeKuttaCoefficients::rungeKuttaFehlberg56)
            .value("rkf_78", tni::RungeKuttaCoefficients::rungeKuttaFehlberg78)
            .value("rkdp_87", tni::RungeKuttaCoefficients::rungeKutta87DormandPrince)
            .export_values();

    py::enum_<tni::ExtrapolationMethodStepSequences>(m, "ExtrapolationMethodStepSequences")
            .value("bulirsch_stoer_sequence", tni::ExtrapolationMethodStepSequences::bulirsch_stoer_sequence)
            .value("deufelhard_sequence", tni::ExtrapolationMethodStepSequences::deufelhard_sequence)
            .export_values();


    py::class_<
            tni::IntegratorSettings<double>,
            std::shared_ptr<tni::IntegratorSettings<double>>>(m, "IntegratorSettings")
            .def(py::init<
                 const tni::AvailableIntegrators,
                 const double,
                 const double,
                 const int,
                 const bool>(),
                 py::arg("integrator_type"),
                 py::arg("initial_time"),
                 py::arg("initial_time_step"),
                 py::arg("save_frequency") = 1,
                 // TODO: Discuss length of this argument: assess_propagation_termination_condition_during_integration_substeps.
                 py::arg("assess_propagation_termination_condition_during_integration_substeps") = false)
            .def_readwrite("initial_time", &tni::IntegratorSettings<double>::initialTime_ );

    m.def("euler",
          &tni::eulerSettings< double >,
          py::arg("initial_time"),
          py::arg("initial_time_step"),
          py::arg("save_frequency") = 1,
          py::arg("assess_termination_on_minor_steps") = false);

    m.def("runge_kutta_4",
          &tni::rungeKutta4Settings< double >,
          py::arg("initial_time"),
          py::arg("initial_time_step"),
          py::arg("save_frequency") = 1,
          py::arg("assess_termination_on_minor_steps") = false);

    m.def("runge_kutta_variable_step_size",
          &tni::rungeKuttaVariableStepSettingsScalarTolerances< double >,
          py::arg("initial_time"),
          py::arg("initial_time_step"),
          py::arg("coefficient_set"),
          py::arg("minimum_step_size"),
          py::arg("maximum_step_size"),
          py::arg("relative_error_tolerance"),
          py::arg("absolute_error_tolerance"),
          py::arg("save_frequency") = 1,
          py::arg("assess_termination_on_minor_steps") = false,
          py::arg("safety_factor") = 0.8,
          py::arg("maximum_factor_increase") = 4.0,
          py::arg("minimum_factor_increase") = 0.1 );

    m.def("runge_kutta_variable_step_size",
          &tni::rungeKuttaVariableStepSettingsVectorTolerances< double >,
          py::arg("initial_time"),
          py::arg("initial_time_step"),
          py::arg("coefficient_set"),
          py::arg("minimum_step_size"),
          py::arg("maximum_step_size"),
          py::arg("relative_error_tolerance"),
          py::arg("absolute_error_tolerance"),
          py::arg("save_frequency") = 1,
          py::arg("assess_termination_on_minor_steps") = false,
          py::arg("safety_factor") = 0.8,
          py::arg("maximum_factor_increase") = 4.0,
          py::arg("minimum_factor_increase") = 0.1);

    m.def("bulirsch_stoer",
          &tni::bulirschStoerIntegratorSettings< double >,
          py::arg("initial_time"),
          py::arg("initial_time_step"),
          py::arg("extrapolation_sequence"),
          py::arg("maximum_number_of_steps"),
          py::arg("minimum_step_size"),
          py::arg("maximum_step_size"),
          py::arg("relative_error_tolerance") = 1.0E-12,
          py::arg("absolute_error_tolerance") = 1.0E-12,
          py::arg("save_frequency") = 1,
          py::arg("check_termination_on_minor_steps") = 0,
          py::arg("safety_factor") = 0.7,
          py::arg("maximum_factor_increase") = 10.0,
          py::arg("minimum_factor_increase") = 0.1 );

    m.def("adams_bashforth_moulton",
          &tni::adamsBashforthMoultonSettings< double >,
          py::arg("initial_time"),
          py::arg("initial_time_step"),
          py::arg("minimum_step_size"),
          py::arg("maximum_step_size"),
          py::arg("relative_error_tolerance") = 1.0E-12,
          py::arg("absolute_error_tolerance") = 1.0E-12,
          py::arg("minimum_order") = 6,
          py::arg("maximum_order") = 11,
          py::arg("save_frequency") = 1,
          py::arg("assess_termination_on_minor_steps") = false,
          py::arg("bandwidth") = 200.0 );
}

void expose_propagator_setup(py::module &m)
{

    py::enum_<tp::TranslationalPropagatorType>(m, "TranslationalPropagatorType")
            .value("undefined_translational_propagator",
                   tp::TranslationalPropagatorType::undefined_translational_propagator)
            .value("cowell",
                   tp::TranslationalPropagatorType::cowell)
            .value("encke",
                   tp::TranslationalPropagatorType::encke)
            .value("gauss_keplerian",
                   tp::TranslationalPropagatorType::gauss_keplerian)
            .value("gauss_modified_equinoctial",
                   tp::TranslationalPropagatorType::gauss_modified_equinoctial)
            .value("unified_state_model_quaternions",
                   tp::TranslationalPropagatorType::unified_state_model_quaternions)
            .value("unified_state_model_modified_rodrigues_parameters",
                   tp::TranslationalPropagatorType::unified_state_model_modified_rodrigues_parameters)
            .value("unified_state_model_exponential_map",
                   tp::unified_state_model_exponential_map)
            .export_values();

    py::enum_<tp::RotationalPropagatorType>(m, "RotationalPropagatorType")
            .value("undefined_rotational_propagator",
                   tp::RotationalPropagatorType::undefined_rotational_propagator)
            .value("quaternions",
                   tp::RotationalPropagatorType::quaternions)
            .value("modified_rodrigues_parameters",
                   tp::RotationalPropagatorType::modified_rodrigues_parameters)
            .value("exponential_map",
                   tp::RotationalPropagatorType::exponential_map)
            .export_values();

    py::class_<tp::DependentVariableSaveSettings,
            std::shared_ptr<tp::DependentVariableSaveSettings>>(m, "DependentVariableSaveSettings")
            .def(py::init<
                 const std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings>>,
                 const bool>(),
                 py::arg("dependent_variables"),
                 py::arg("print_dependent_variable_types") = true);

    py::class_<
            tp::PropagatorSettings<double>,
            std::shared_ptr<tp::PropagatorSettings<double>>>(m, "PropagatorSettings")
            .def("reset_initial_states", &tp::PropagatorSettings<double>::resetInitialStates);

    py::class_<
            tp::SingleArcPropagatorSettings<double>,
            std::shared_ptr<tp::SingleArcPropagatorSettings<double>>,
            tp::PropagatorSettings<double>>(m, "SingleArcPropagatorSettings")
            .def_property("termination_settings",
                          &tp::SingleArcPropagatorSettings<double>::getTerminationSettings,
                          &tp::SingleArcPropagatorSettings<double>::resetTerminationSettings);
    py::class_<
            tp::TranslationalStatePropagatorSettings<double>,
            std::shared_ptr<tp::TranslationalStatePropagatorSettings<double>>,
            tp::SingleArcPropagatorSettings<double>>(m, "TranslationalStatePropagatorSettings")
            .def(// ctor 1
                 py::init<
                 const std::vector<std::string> &,
                 const tba::AccelerationMap &,
                 const std::vector<std::string> &,
                 const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
                 const std::shared_ptr<tp::PropagationTerminationSettings>,
                 const tp::TranslationalPropagatorType,
                 const std::shared_ptr<tp::DependentVariableSaveSettings>,
                 const double>(),
                 py::arg("central_bodies"),
                 py::arg("acceleration_models"),
                 py::arg("bodies_to_integrate"),
                 py::arg("initial_states"),
                 py::arg("termination_settings"),
                 py::arg("propagator") = tp::TranslationalPropagatorType::cowell,
                 py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
                 py::arg("print_interval") = TUDAT_NAN)
            .def(// ctor 2
                 py::init<const std::vector<std::string> &,
                 const tss::SelectedAccelerationMap &,
                 const std::vector<std::string> &,
                 const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
                 const std::shared_ptr<tp::PropagationTerminationSettings>,
                 const tp::TranslationalPropagatorType,
                 const std::shared_ptr<tp::DependentVariableSaveSettings>,
                 const double>(),
                 py::arg("central_bodies"),
                 py::arg("acceleration_settings"),
                 py::arg("bodies_to_integrate"),
                 py::arg("initial_states"),
                 py::arg("termination_settings"),
                 py::arg("propagator") = tp::cowell,
                 py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
                 py::arg("print_interval") = TUDAT_NAN)
            .def(// ctor 3
                 py::init<const std::vector<std::string> &,
                 const tba::AccelerationMap &,
                 const std::vector<std::string> &,
                 const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
                 const double,
                 const tp::TranslationalPropagatorType,
                 const std::shared_ptr<tp::DependentVariableSaveSettings>,
                 const double>(),
                 py::arg("central_bodies"),
                 py::arg("acceleration_models"),
                 py::arg("bodies_to_integrate"),
                 py::arg("initial_states"),
                 py::arg("termination_time"),
                 py::arg("propagator") = tp::cowell,
                 py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
                 py::arg("print_interval") = TUDAT_NAN)
            .def(// ctor 4
                 py::init<const std::vector<std::string> &,
                 const tss::SelectedAccelerationMap &,
                 const std::vector<std::string> &,
                 const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
                 const double,
                 const tp::TranslationalPropagatorType,
                 const std::shared_ptr<tp::DependentVariableSaveSettings>,
                 const double>(),
                 py::arg("central_bodies"),
                 py::arg("acceleration_settings"),
                 py::arg("bodies_to_integrate"),
                 py::arg("initial_states"),
                 py::arg("termination_time"),
                 py::arg("propagator") = tp::cowell,
                 py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
                 py::arg("print_interval") = TUDAT_NAN)
            .def("recreate_state_derivative_models", &tp::TranslationalStatePropagatorSettings<double>::resetIntegratedStateModels,
                 py::arg("bodies") )
            .def("get_propagated_state_size", &tp::TranslationalStatePropagatorSettings<double>::getPropagatedStateSize)
            .def("reset_and_recreate_acceleration_models", &tp::TranslationalStatePropagatorSettings<double>::resetAccelerationModelsMap,
                 py::arg("new_acceleration_settings"),
                 py::arg("bodies") )
            .def_property_readonly("acceleration_settings", &tp::TranslationalStatePropagatorSettings<double>::getAccelerationSettingsMap);


    py::class_<
            tp::MultiTypePropagatorSettings<double>,
            std::shared_ptr<tp::MultiTypePropagatorSettings<double>>,
            tp::SingleArcPropagatorSettings<double>>(m, "MultiTypePropagatorSettings")
            .def("reset_initial_states", &tp::MultiTypePropagatorSettings<double>::resetInitialStates)
            .def("recreate_state_derivative_models", &tp::MultiTypePropagatorSettings<double>::resetIntegratedStateModels,
                 py::arg("bodies") )
            .def("single_type_settings", &tp::MultiTypePropagatorSettings<double>::getSingleTypePropagatorSettings,
                 py::arg("state_type") )
            .def_property_readonly("propagator_settings_per_type", &tp::MultiTypePropagatorSettings<double>::getPropagatorSettingsMap);




    py::class_<
            tp::MassPropagatorSettings<double>,
            std::shared_ptr<tp::MassPropagatorSettings<double>>,
            tp::SingleArcPropagatorSettings<double>>(m, "MassPropagatorSettings");


    m.def("combine_initial_states",
          &tp::createCombinedInitialState<double>,
          py::arg("propagator_settings_per_type") );

    m.def("translational",
          py::overload_cast<
          const std::vector< std::string >&,
          const tba::AccelerationMap&,
          const std::vector< std::string >&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const std::shared_ptr<tp::PropagationTerminationSettings>,
          const tp::TranslationalPropagatorType,
          const std::shared_ptr<tp::DependentVariableSaveSettings>,
          const double>(&tp::translationalStatePropagatorSettings<double>),
          py::arg("central_bodies"),
          py::arg("acceleration_models"),
          py::arg("bodies_to_integrate"),
          py::arg("initial_states"),
          py::arg("termination_settings"),
          py::arg("propagator") = tp::cowell,
          py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
          py::arg("print_interval") = TUDAT_NAN);


    m.def("translational",
          py::overload_cast<
          const std::vector< std::string >&,
          const tba::AccelerationMap&,
          const std::vector< std::string >&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const std::shared_ptr<tp::PropagationTerminationSettings>,
          const tp::TranslationalPropagatorType,
          const std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >&,
          const double>(&tp::translationalStatePropagatorSettings<double>),
          py::arg("central_bodies"),
          py::arg("acceleration_models"),
          py::arg("bodies_to_integrate"),
          py::arg("initial_states"),
          py::arg("termination_settings"),
          py::arg("propagator") = tp::cowell,
          py::arg("output_variables") =  std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >(),
          py::arg("print_interval") = TUDAT_NAN);


    m.def("translational",
          py::overload_cast<
          const std::vector< std::string >&,
          const tba::AccelerationMap&,
          const std::vector< std::string >&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const double,
          const tp::TranslationalPropagatorType,
          const std::shared_ptr<tp::DependentVariableSaveSettings>,
          const double>(&tp::translationalStatePropagatorSettings<double>),
          py::arg("central_bodies"),
          py::arg("acceleration_models"),
          py::arg("bodies_to_integrate"),
          py::arg("initial_states"),
          py::arg("termination_time"),
          py::arg("propagator") = tp::cowell,
          py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
          py::arg("print_interval") = TUDAT_NAN);

    m.def("translational",
          py::overload_cast<
          const std::vector< std::string >&,
          const tba::AccelerationMap&,
          const std::vector< std::string >&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const double,
          const tp::TranslationalPropagatorType,
          const std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >&,
          const double>(&tp::translationalStatePropagatorSettings<double>),
          py::arg("central_bodies"),
          py::arg("acceleration_models"),
          py::arg("bodies_to_integrate"),
          py::arg("initial_states"),
          py::arg("termination_time"),
          py::arg("propagator") = tp::cowell,
          py::arg("output_variables") =  std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >(),
          py::arg("print_interval") = TUDAT_NAN);

    m.def("translational",
          py::overload_cast<
          const std::vector< std::string >&,
          const tss::SelectedAccelerationMap&,
          const std::vector< std::string >&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const double,
          const tp::TranslationalPropagatorType,
          const std::shared_ptr<tp::DependentVariableSaveSettings>,
          const double>(&tp::translationalStatePropagatorSettings<double>),
          py::arg("central_bodies"),
          py::arg("acceleration_settings"),
          py::arg("bodies_to_integrate"),
          py::arg("initial_states"),
          py::arg("termination_time"),
          py::arg("propagator") = tp::cowell,
          py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
          py::arg("print_interval") = TUDAT_NAN);

    m.def("translational",
          py::overload_cast<
          const std::vector< std::string >&,
          const tss::SelectedAccelerationMap&,
          const std::vector< std::string >&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const std::shared_ptr<tp::PropagationTerminationSettings>,
          const tp::TranslationalPropagatorType,
          const std::shared_ptr<tp::DependentVariableSaveSettings>,
          const double>(&tp::translationalStatePropagatorSettings<double>),
          py::arg("central_bodies"),
          py::arg("acceleration_settings"),
          py::arg("bodies_to_integrate"),
          py::arg("initial_states"),
          py::arg("termination_settings"),
          py::arg("propagator") = tp::cowell,
          py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
          py::arg("print_interval") = TUDAT_NAN);

    m.def("translational",
          py::overload_cast<
          const std::vector< std::string >&,
          const tss::SelectedAccelerationMap&,
          const std::vector< std::string >&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const std::shared_ptr<tp::PropagationTerminationSettings>,
          const tp::TranslationalPropagatorType,
          const std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >&,
          const double>(&tp::translationalStatePropagatorSettings<double>),
          py::arg("central_bodies"),
          py::arg("acceleration_settings"),
          py::arg("bodies_to_integrate"),
          py::arg("initial_states"),
          py::arg("termination_settings"),
          py::arg("propagator") = tp::cowell,
          py::arg("output_variables") = std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >(),
          py::arg("print_interval") = TUDAT_NAN);

    m.def("mass",
          py::overload_cast<
          const std::vector< std::string >,
          const std::map< std::string, std::shared_ptr< tba::MassRateModel > >&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const std::shared_ptr< tp::PropagationTerminationSettings >,
          const std::shared_ptr< tp::DependentVariableSaveSettings >,
          const double >(&tp::massPropagatorSettings<double>),
          py::arg("bodies_with_mass_to_propagate"),
          py::arg("mass_rate_models"),
          py::arg("initial_body_masses"),
          py::arg("termination_settings"),
          py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
          py::arg("print_interval") = TUDAT_NAN);

    m.def("mass",
          py::overload_cast<
          const std::vector< std::string >,
          const std::map< std::string, std::vector< std::shared_ptr< tba::MassRateModel > > >&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const std::shared_ptr< tp::PropagationTerminationSettings >,
          const std::shared_ptr< tp::DependentVariableSaveSettings >,
          const double >(&tp::massPropagatorSettings<double>),
          py::arg("bodies_with_mass_to_propagate"),
          py::arg("mass_rate_models"),
          py::arg("initial_body_masses"),
          py::arg("termination_settings"),
          py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
          py::arg("print_interval") = TUDAT_NAN);

    m.def("rotational",
          py::overload_cast<
          const tba::TorqueModelMap&,
          const std::vector< std::string >&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const std::shared_ptr< tp::PropagationTerminationSettings >,
          const tp::RotationalPropagatorType,
          const std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >,
          const double
          >(&tp::rotatonalPropagatorSettings<double>),
          py::arg("torque_models"),
          py::arg("bodies_to_propagate"),
          py::arg("initial_states"),
          py::arg("termination_settings"),
          py::arg("propagator") = tp::quaternions,
          py::arg("output_variables") = std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >(),
          py::arg("print_interval") = TUDAT_NAN);


    m.def("rotational",
          py::overload_cast<
          const tss::SelectedTorqueMap&,
          const std::vector< std::string >&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const std::shared_ptr< tp::PropagationTerminationSettings >,
          const tp::RotationalPropagatorType,
          const std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >,
          const double
          >(&tp::rotatonalPropagatorSettings<double>),
          py::arg("torque_model_settings"),
          py::arg("bodies_to_propagate"),
          py::arg("initial_states"),
          py::arg("termination_settings"),
          py::arg("propagator") = tp::quaternions,
          py::arg("output_variables") = std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >(),
          py::arg("print_interval") = TUDAT_NAN);


    m.def("mass",
          py::overload_cast<
          const std::vector< std::string >,
          const tss::SelectedMassRateModelMap&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const std::shared_ptr< tp::PropagationTerminationSettings >,
          const std::shared_ptr< tp::DependentVariableSaveSettings >,
          const double >(&tp::massPropagatorSettings<double>),
          py::arg("bodies_with_mass_to_propagate"),
          py::arg("mass_rate_settings"),
          py::arg("initial_body_masses"),
          py::arg("termination_settings"),
          py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
          py::arg("print_interval") = TUDAT_NAN);

    m.def("mass",
          py::overload_cast<
          const std::vector< std::string >,
          const std::map< std::string, std::shared_ptr< tba::MassRateModel > >&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const std::shared_ptr< tp::PropagationTerminationSettings >,
          const std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >,
          const double >(&tp::massPropagatorSettings<double>),
          py::arg("bodies_with_mass_to_propagate"),
          py::arg("mass_rate_models"),
          py::arg("initial_body_masses"),
          py::arg("termination_settings"),
          py::arg("output_variables") = std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >(),
          py::arg("print_interval") = TUDAT_NAN);

    m.def("mass",
          py::overload_cast<
          const std::vector< std::string >,
          const std::map< std::string, std::vector< std::shared_ptr< tba::MassRateModel > > >&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const std::shared_ptr< tp::PropagationTerminationSettings >,
          const std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >,
          const double >(&tp::massPropagatorSettings<double>),
          py::arg("bodies_with_mass_to_propagate"),
          py::arg("mass_rate_models"),
          py::arg("initial_body_masses"),
          py::arg("termination_settings"),
          py::arg("output_variables") = std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >(),
          py::arg("print_interval") = TUDAT_NAN);

    m.def("mass",
          py::overload_cast<
          const std::vector< std::string >,
          const tss::SelectedMassRateModelMap&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const std::shared_ptr< tp::PropagationTerminationSettings >,
          const std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >,
          const double >(&tp::massPropagatorSettings<double>),
          py::arg("bodies_with_mass_to_propagate"),
          py::arg("mass_rate_settings"),
          py::arg("initial_body_masses"),
          py::arg("termination_settings"),
          py::arg("output_variables") = std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >(),
          py::arg("print_interval") = TUDAT_NAN);


    m.def("multitype",
          py::overload_cast<
          const std::vector< std::shared_ptr< tp::SingleArcPropagatorSettings< double > > >,
          const std::shared_ptr< tp::PropagationTerminationSettings >,
          const std::shared_ptr< tp::DependentVariableSaveSettings >,
          const double >( &tp::multiTypePropagatorSettings<double> ),
          py::arg("propagator_settings_list"),
          py::arg("termination_settings"),
          py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>( ),
          py::arg("print_interval") = TUDAT_NAN );


    m.def("multitype",
          py::overload_cast<
          const std::vector< std::shared_ptr< tp::SingleArcPropagatorSettings< double > > >,
          const std::shared_ptr< tp::PropagationTerminationSettings >,
          const std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >,
          const double >( &tp::multiTypePropagatorSettings<double> ),
          py::arg("propagator_settings_list"),
          py::arg("termination_settings"),
          py::arg("output_variables") = std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >( ),
          py::arg("print_interval") = TUDAT_NAN );

    m.def("multi_arc",
          &tp::multiArcPropagatorSettings<double>,
          py::arg("single_arc_settings"),
          py::arg("transfer_state_to_next_arc") = false );

    m.def("hybrid_arc",
          &tp::hybridArcPropagatorSettings<double>,
          py::arg("single_arc_settings"),
          py::arg("multi_arc_settings") );

    py::class_<tp::PropagationTerminationSettings,
            std::shared_ptr<tp::PropagationTerminationSettings>>
            PropagationTerminationSettings_(m, "PropagationTerminationSettings");

    py::class_<
            tp::PropagationDependentVariableTerminationSettings,
            std::shared_ptr<tp::PropagationDependentVariableTerminationSettings>,
            tp::PropagationTerminationSettings>(m, "PropagationDependentVariableTerminationSettings")
            .def(py::init<
                 const std::shared_ptr<tp::SingleDependentVariableSaveSettings>,
                 const double,
                 const bool,
                 const bool,
                 const std::shared_ptr<tudat::root_finders::RootFinderSettings>>(),
                 py::arg("dependent_variadble_settings"),
                 py::arg("limit_value"),
                 py::arg("use_as_lower_limit"),
                 py::arg("terminate_exactly_on_final_condition") = false,
                 py::arg("termination_root_finder_settings") = nullptr);

    m.def("time_termination",
          &tp::propagationTimeTerminationSettings,
          py::arg("termination_time"),
          py::arg("terminate_exactly_on_final_condition") = false);

    m.def("cpu_time_termination",
          &tp::propagationCPUTimeTerminationSettings,
          py::arg("cpu_termination_time") );

    m.def("dependent_variable_termination",
          &tp::propagationDependentVariableTerminationSettings,
          py::arg("dependent_variable_settings"),
          py::arg("limit_value"),
          py::arg("use_as_lower_limit"),
          py::arg("terminate_exactly_on_final_condition") = false,
          py::arg("termination_root_finder_settings") = nullptr );

    m.def("custom_termination",
          &tp::popagationCustomTerminationSettings,
          py::arg("custom_condition"));


    m.def("hybrid_termination",
          &tp::propagationHybridTerminationSettings,
          py::arg("termination_settings"),
          py::arg("fulfill_single_condition") );
}

void expose_propagation_setup(py::module &m) {

    py::enum_<tss::ThrustFrames>(m, "ThrustFrames")
            .value("unspecified_thrust_frame", tss::ThrustFrames::unspecified_thrust_frame)
            .value("inertial_thurst_frame", tss::ThrustFrames::inertial_thurst_frame)
            .value("lvlh_thrust_frame", tss::ThrustFrames::lvlh_thrust_frame)
            .export_values();
    /*
   * propagation_setup
   *  ├── accelerationSettings.h
   *  ├── createAccelerationModels.h
   *  ├── createEnvironmentUpdater.h
   *  ├── createMassRateModels.h
   *  ├── createStateDerivativeModel.h
   *  ├── createThrustModelGuidance.h
   *  ├── createTorqueModel.h
   *  ├── dynamicsSimulator.h
   *  ├── environmentUpdater.h
   *  ├── propagationCR3BPFullProblem.h
   *  ├── propagationLambertTargeterFullProblem.h
   *  ├── propagationOutput.h
   *  ├── propagationOutputSettings.h
   *  ├── propagationPatchedConicFullProblem.h
   *  ├── propagationSettings.h
   *  ├── propagationTermination.h
   *  ├── propagationTerminationSettings.h
   *  ├── setNumericallyIntegratedStates.h
   *  ├── thrustSettings.h
   *  └── torqueSettings.h
   *
   * propagation_setup/
   *  ├── createAccelerationModels.cpp
   *  ├── createEnvironmentUpdater.cpp
   *  ├── createMassRateModels.cpp
   *  ├── createStateDerivativeModel.cpp
   *  ├── createThrustModelGuidance.cpp
   *  ├── createTorqueModel.cpp
   *  ├── dynamicsSimulator.cpp
   *  ├── environmentUpdater.cpp
   *  ├── propagationCR3BPFullProblem.cpp
   *  ├── propagationLambertTargeterFullProblem.cpp
   *  ├── propagationOutput.cpp
   *  ├── propagationOutputSettings.cpp
   *  ├── propagationPatchedConicFullProblem.cpp
   *  ├── propagationSettings.cpp
   *  ├── propagationTermination.cpp
   *  ├── setNumericallyIntegratedStates.cpp
   *  └── thrustSettings.cpp
   *
   */


    //////////////////////////////////////////////////////////////////////////////
    // propagationTerminationSettings.h
    //////////////////////////////////////////////////////////////////////////////
    //  py::enum_<tss::PropagationTerminationTypes,
    //            std::shared_ptr<>>
    //  enum PropagationTerminationTypes
    //  {
    //    time_stopping_condition = 0,
    //    cpu_time_stopping_condition = 1,
    //    dependent_variable_stopping_condition = 2,
    //    hybrid_stopping_condition = 3,
    //    custom_stopping_condition = 4
    //  };

    //  py::class_<tss::ThrustAccelerationSettings,
    //             std::shared_ptr<tss::ThrustAccelerationSettings>,
    //             tss::AccelerationSettings>(m, "ThrustAccelerationSettings")
    //      .def(py::init<//ctor 1
    //               const std::shared_ptr<tss::ThrustDirectionGuidanceSettings>,
    //               const std::shared_ptr<tss::ThrustMagnitudeSettings>>(),
    //           py::arg("thrust_direction_settings"),
    //           py::arg("thrust_magnitude_settings"));

    //////////////////////////////////////////////////////////////////////////////
    // createAccelerationModels.cpp
    //////////////////////////////////////////////////////////////////////////////
    m.def("create_acceleration_models",// overload [1/2]
          py::overload_cast<const tss::SystemOfBodies &,
          const tss::SelectedAccelerationMap &,
          const std::vector<std::string> &,
          const std::vector<std::string> &>(
              &tss::createAccelerationModelsMap),
          py::arg("body_system"),
          py::arg("selected_acceleration_per_body"),
          py::arg("bodies_to_propagate"),
          py::arg("central_bodies"));

    m.def("create_acceleration_models",// overload [2/2]
          py::overload_cast<const tss::SystemOfBodies &,
          const tss::SelectedAccelerationMap &,
          const std::map<std::string, std::string> &>(
              &tss::createAccelerationModelsMap),
          py::arg("body_system"),
          py::arg("selected_acceleration_per_body"),
          py::arg("central_bodies"));

    //////////////////////////////////////////////////////////////////////////////
    // createTorqueModels.cpp
    //////////////////////////////////////////////////////////////////////////////

    m.def("create_torque_models",// overload [1/2]
              &tss::createTorqueModelsMap,
          py::arg("body_system"),
          py::arg("selected_acceleration_per_body"),
          py::arg("bodies_to_propagate"));


    //////////////////////////////////////////////////////////////////////////////
    // dynamicsSimulator.h / dynamicsSimulator.cpp
    //////////////////////////////////////////////////////////////////////////////
    //  m.def("get_initial_state_of_bodies",// overload [1/2]
    //        py::overload_cast<const std::vector<std::string> &,
    //                          const std::vector<std::string> &,
    //                          const tss::SystemOfBodies &,
    //                          const double,
    //                          std::shared_ptr<te::ReferenceFrameManager>>(
    //            &tp::getInitialStatesOfBodies<>));

    py::enum_<tp::IntegratedStateType>(m, "StateType")
            .value("hybrid_type", tp::IntegratedStateType::hybrid)
            .value("translational_type", tp::IntegratedStateType::translational_state)
            .value("rotational_type", tp::IntegratedStateType::rotational_state)
            .value("mass_type", tp::IntegratedStateType::body_mass_state)
            .value("custom_type", tp::IntegratedStateType::custom_state)
            .export_values();

    m.def("get_initial_state_of_bodies",// overload [2/2]
          py::overload_cast<const std::vector<std::string> &,
          const std::vector<std::string> &,
          const tss::SystemOfBodies &,
          const double>(
              &tp::getInitialStatesOfBodies<>),
          py::arg("bodies_to_propagate"),
          py::arg("central_bodies"),
          py::arg("body_system"),
          py::arg("initial_time"));

    py::class_<
            tp::SingleArcDynamicsSimulator<double, double>,
            std::shared_ptr<tp::SingleArcDynamicsSimulator<double, double>>>(m, "SingleArcDynamicsSimulator")
            .def(py::init<
                 const tudat::simulation_setup::SystemOfBodies &,
                 const std::shared_ptr<tudat::numerical_integrators::IntegratorSettings<double>>,
                 const std::shared_ptr<tp::PropagatorSettings<double>>,
                 const bool,
                 const bool,
                 const bool,
                 const bool,
                 const std::chrono::steady_clock::time_point,
                 const std::vector<std::shared_ptr<tp::SingleStateTypeDerivative<double, double>>>&,
                 const bool >(),
                 py::arg("body_map"),
                 py::arg("integrator_settings"),
                 py::arg("propagator_settings"),
                 py::arg("are_equations_of_motion_to_be_integrated") = true,
                 py::arg("clear_numerical_solutions") = false,
                 py::arg("set_integrated_result") = false,
                 py::arg("print_number_of_function_evaluations") = false,
                 py::arg("initial_clock_time") = std::chrono::steady_clock::now(),
                 py::arg("state_derivative_models") =
            std::vector<std::shared_ptr<tp::SingleStateTypeDerivative<double, double>>>(),
                 py::arg("print_dependent_variable_data" )= true )
            .def("integrate_equations_of_motion",
                 &tp::SingleArcDynamicsSimulator<double, double>::integrateEquationsOfMotion,
                 py::arg("initial_states"))
            .def("get_equations_of_motion_numerical_solution",
                 &tp::SingleArcDynamicsSimulator<double, double>::getEquationsOfMotionNumericalSolution)
            .def_property_readonly("state_history",
                                   &tp::SingleArcDynamicsSimulator<double, double>::getEquationsOfMotionNumericalSolution)
            .def("get_equations_of_motion_numerical_solution_raw",
                 &tp::SingleArcDynamicsSimulator<double, double>::getEquationsOfMotionNumericalSolutionRaw)
            .def_property_readonly("unprocessed_state_history",
                                   &tp::SingleArcDynamicsSimulator<double, double>::getEquationsOfMotionNumericalSolutionRaw)
            .def("get_dependent_variable_history",
                 &tp::SingleArcDynamicsSimulator<double, double>::getDependentVariableHistory)
            .def_property_readonly("dependent_variable_history",
                                   &tp::SingleArcDynamicsSimulator<double, double>::getDependentVariableHistory)
            .def("get_cumulative_computation_time_history",
                 &tp::SingleArcDynamicsSimulator<double, double>::getCumulativeComputationTimeHistory)
            .def("get_cumulative_number_of_function_evaluations",
                 &tp::SingleArcDynamicsSimulator<double, double>::getCumulativeNumberOfFunctionEvaluations)
            .def("get_equations_of_motion_numerical_solution_base",
                 &tp::SingleArcDynamicsSimulator<double, double>::getEquationsOfMotionNumericalSolutionBase)
            .def("get_dependent_variable_numerical_solution_base",
                 &tp::SingleArcDynamicsSimulator<double, double>::getDependentVariableNumericalSolutionBase)
            .def("get_cumulative_computation_time_history_base",
                 &tp::SingleArcDynamicsSimulator<double, double>::getCumulativeComputationTimeHistoryBase)
            .def("manually_set_and_process_raw_numerical_equations_of_motion_solution",
                 &tp::SingleArcDynamicsSimulator<double, double>::manuallySetAndProcessRawNumericalEquationsOfMotionSolution,
                 py::arg("equations_of_motion_numerical_solution"),
                 py::arg("dependent_variable_history"),
                 py::arg("process_solution"))
            .def("get_integrator_settings",
                 &tp::SingleArcDynamicsSimulator<double, double>::getIntegratorSettings)
            .def("get_state_derivative_function",
                 &tp::SingleArcDynamicsSimulator<double, double>::getStateDerivativeFunction)
            .def("get_double_state_derivative_function",
                 &tp::SingleArcDynamicsSimulator<double, double>::getDoubleStateDerivativeFunction)
            .def("get_environment_updater",
                 &tp::SingleArcDynamicsSimulator<double, double>::getEnvironmentUpdater)
            .def("get_dynamics_state_derivative",
                 &tp::SingleArcDynamicsSimulator<double, double>::getDynamicsStateDerivative)
            .def("get_propagation_termination_condition",
                 &tp::SingleArcDynamicsSimulator<double, double>::getPropagationTerminationCondition)
            .def("get_integrated_state_processors",
                 &tp::SingleArcDynamicsSimulator<double, double>::getIntegratedStateProcessors)
            .def("get_propagation_termination_reason",
                 &tp::SingleArcDynamicsSimulator<double, double>::getPropagationTerminationReason)
            .def("integration_completed_successfully",
                 &tp::SingleArcDynamicsSimulator<double, double>::integrationCompletedSuccessfully)
            .def("get_dependent_variable_ids",
                 &tp::SingleArcDynamicsSimulator<double, double>::getDependentVariableIds)
            .def("get_initial_propagation_time",
                 &tp::SingleArcDynamicsSimulator<double, double>::getInitialPropagationTime)
            .def("reset_initial_propagation_time",
                 &tp::SingleArcDynamicsSimulator<double, double>::resetInitialPropagationTime)
            .def("get_dependent_variables_functions",
                 &tp::SingleArcDynamicsSimulator<double, double>::getDependentVariablesFunctions)
            .def("reset_propagation_termination_conditions",
                 &tp::SingleArcDynamicsSimulator<double, double>::resetPropagationTerminationConditions)
            .def("process_numerical_equations_of_motion_solution",
                 &tp::SingleArcDynamicsSimulator<double, double>::processNumericalEquationsOfMotionSolution)
            .def("suppress_dependent_variable_terminal_printing",
                 &tp::SingleArcDynamicsSimulator<double, double>::suppressDependentVariableDataPrinting)
            .def("enable_dependent_variable_terminal_printing",
                 &tp::SingleArcDynamicsSimulator<double, double>::enableDependentVariableDataPrinting);


    //        py::enum_<tp::VariableType>(m, "VariableType")
    //                .value("independent_variable", tp::VariableType::independentVariable)
    //                .value("cpu_time_variable", tp::VariableType::cpuTimeVariable)
    //                .value("state_variable", tp::VariableType::stateVariable)
    //                .value("dependent_variable", tp::VariableType::dependentVariable)
    //                .export_values();



    auto acceleration_setup = m.def_submodule("acceleration");
    expose_acceleration_setup(acceleration_setup);

    auto torque_setup = m.def_submodule("torque");
    expose_torque_setup(torque_setup);

    auto integrator_setup = m.def_submodule("integrator");
    expose_integrator_setup(integrator_setup);

    auto propagator_setup = m.def_submodule("propagator");
    expose_propagator_setup(propagator_setup);

    auto mass_setup = m.def_submodule("mass");
    expose_mass_rate_setup(mass_setup);


    auto dependent_variable_setup = m.def_submodule("dependent_variable");
    expose_dependent_variable_setup(dependent_variable_setup);
}

}// namespace tudatpy
