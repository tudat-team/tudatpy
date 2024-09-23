/*    Copyright (c) 2010-2021, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define PYBIND11_DETAILED_ERROR_MESSAGES

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/simulation/estimation_setup/createObservationModel.h>
#include <tudat/simulation/estimation_setup/observationSimulationSettings.h>
#include <tudat/simulation/estimation_setup/processOdfFile.h>
#include <tudat/simulation/estimation_setup/processTrackingTxtFile.h>
#include <tudat/simulation/estimation_setup/simulateObservations.h>

#include "tudatpy/scalarTypes.h"


namespace py = pybind11;

namespace tss = tudat::simulation_setup;
namespace tom = tudat::observation_models;
namespace tuc = tudat::unit_conversions;

namespace tudat {

    namespace simulation_setup {


        void addNoiseFunctionToObservationSimulationSettingsPy(
            const std::vector<
                std::shared_ptr<ObservationSimulationSettings<TIME_TYPE>>>&
                observationSimulationSettings,
            const std::function<Eigen::VectorXd(const double)>
                observationNoiseFunction) {
            tss::addNoiseFunctionToObservationSimulationSettings<
                TIME_TYPE, Eigen::VectorXd>(observationSimulationSettings,
                                            observationNoiseFunction);
        }


        void addNoiseFunctionToObservationSimulationSettingsPy(
            const std::vector<
                std::shared_ptr<ObservationSimulationSettings<TIME_TYPE>>>&
                observationSimulationSettings,
            const std::function<Eigen::VectorXd(const double)>
                observationNoiseFunction,
            const tom::ObservableType observableType) {
            tss::addNoiseFunctionToObservationSimulationSettings<
                TIME_TYPE, Eigen::VectorXd, const tom::ObservableType>(
                observationSimulationSettings, observationNoiseFunction,
                observableType);
        }


        void addNoiseFunctionToObservationSimulationSettingsPy(
            const std::vector<
                std::shared_ptr<ObservationSimulationSettings<TIME_TYPE>>>&
                observationSimulationSettings,
            const std::function<Eigen::VectorXd(const double)>
                observationNoiseFunction,
            const tom::ObservableType observableType,
            const tom::LinkDefinition& linkEnds) {
            tss::addNoiseFunctionToObservationSimulationSettings<
                TIME_TYPE, Eigen::VectorXd, const tom::ObservableType,
                const tom::LinkDefinition&>(observationSimulationSettings,
                                            observationNoiseFunction,
                                            observableType, linkEnds);
        }

        void addGaussianNoiseFunctionToObservationSimulationSettingsPy(
            const std::vector<
                std::shared_ptr<ObservationSimulationSettings<TIME_TYPE>>>&
                observationSimulationSettings,
            const double observationNoiseAmplitude) {
            tss::addGaussianNoiseFunctionToObservationSimulationSettings<
                TIME_TYPE>(observationSimulationSettings,
                           observationNoiseAmplitude);
        }


        void addGaussianNoiseFunctionToObservationSimulationSettingsPy(
            const std::vector<
                std::shared_ptr<ObservationSimulationSettings<TIME_TYPE>>>&
                observationSimulationSettings,
            const double observationNoiseAmplitude,
            const tom::ObservableType observableType) {
            tss::addGaussianNoiseFunctionToObservationSimulationSettings<
                TIME_TYPE, const tom::ObservableType>(
                observationSimulationSettings, observationNoiseAmplitude,
                observableType);
        }


        void addGaussianNoiseFunctionToObservationSimulationSettingsPy(
            const std::vector<
                std::shared_ptr<ObservationSimulationSettings<TIME_TYPE>>>&
                observationSimulationSettings,
            const double observationNoiseAmplitude,
            const tom::ObservableType observableType,
            const tom::LinkDefinition& linkEnds) {
            tss::addGaussianNoiseFunctionToObservationSimulationSettings<
                TIME_TYPE, const tom::ObservableType,
                const tom::LinkDefinition&>(observationSimulationSettings,
                                            observationNoiseAmplitude,
                                            observableType, linkEnds);
        }

        void addViabilityToObservationSimulationSettingsPy(
            const std::vector<
                std::shared_ptr<ObservationSimulationSettings<TIME_TYPE>>>&
                observationSimulationSettings,
            const std::vector<std::shared_ptr<
                observation_models::ObservationViabilitySettings>>&
                viabilitySettingsList) {
            tss::addViabilityToObservationSimulationSettings<TIME_TYPE>(
                observationSimulationSettings, viabilitySettingsList);
        }

        void addViabilityToObservationSimulationSettingsPy(
            const std::vector<
                std::shared_ptr<ObservationSimulationSettings<TIME_TYPE>>>&
                observationSimulationSettings,
            const std::vector<std::shared_ptr<
                observation_models::ObservationViabilitySettings>>&
                viabilitySettingsList,
            const tom::ObservableType observableType) {
            tss::addViabilityToObservationSimulationSettings<
                TIME_TYPE, const tom::ObservableType>(
                observationSimulationSettings, viabilitySettingsList,
                observableType);
        }

        void addViabilityToObservationSimulationSettingsPy(
            const std::vector<
                std::shared_ptr<ObservationSimulationSettings<TIME_TYPE>>>&
                observationSimulationSettings,
            const std::vector<std::shared_ptr<
                observation_models::ObservationViabilitySettings>>&
                viabilitySettingsList,
            const tom::ObservableType observableType,
            const tom::LinkDefinition& linkEnds) {
            tss::addViabilityToObservationSimulationSettings<
                TIME_TYPE, const tom::ObservableType,
                const tom::LinkDefinition&>(observationSimulationSettings,
                                            viabilitySettingsList,
                                            observableType, linkEnds);
        }

        void addDependentVariablesToObservationSimulationSettingsPy(
            const std::vector<
                std::shared_ptr<ObservationSimulationSettings<TIME_TYPE>>>&
                observationSimulationSettings,
            const std::vector<
                std::shared_ptr<ObservationDependentVariableSettings>>&
                dependentVariableList,
            const SystemOfBodies& bodies) {
            tss::addDependentVariablesToObservationSimulationSettings<
                TIME_TYPE>(observationSimulationSettings, dependentVariableList,
                           bodies);
        }

        void addDependentVariablesToObservationSimulationSettingsPy(
            const std::vector<
                std::shared_ptr<ObservationSimulationSettings<TIME_TYPE>>>&
                observationSimulationSettings,
            const std::vector<
                std::shared_ptr<ObservationDependentVariableSettings>>&
                dependentVariableList,
            const SystemOfBodies& bodies,
            const tom::ObservableType observableType) {
            tss::addDependentVariablesToObservationSimulationSettings<
                TIME_TYPE, const tom::ObservableType>(
                observationSimulationSettings, dependentVariableList, bodies,
                observableType);
        }

        void addDependentVariablesToObservationSimulationSettingsPy(
            const std::vector<
                std::shared_ptr<ObservationSimulationSettings<TIME_TYPE>>>&
                observationSimulationSettings,
            const std::vector<
                std::shared_ptr<ObservationDependentVariableSettings>>&
                dependentVariableList,
            const SystemOfBodies& bodies,
            const tom::ObservableType observableType,
            const tom::LinkDefinition& linkEnds) {
            tss::addDependentVariablesToObservationSimulationSettings<
                TIME_TYPE, const tom::ObservableType,
                const tom::LinkDefinition&>(observationSimulationSettings,
                                            dependentVariableList, bodies,
                                            observableType, linkEnds);
        }

        void addAncilliarySettingsToObservationSimulationSettingsPy(
            const std::vector<
                std::shared_ptr<ObservationSimulationSettings<TIME_TYPE>>>&
                observationSimulationSettings,
            const std::shared_ptr<tom::ObservationAncilliarySimulationSettings>&
                ancilliarySettings,
            const tom::ObservableType observableType) {
            tss::addAncilliarySettingsToObservationSimulationSettings<
                TIME_TYPE, const tom::ObservableType>(
                observationSimulationSettings, ancilliarySettings,
                observableType);
        }

        void addAncilliarySettingsToObservationSimulationSettingsPy(
            const std::vector<
                std::shared_ptr<ObservationSimulationSettings<TIME_TYPE>>>&
                observationSimulationSettings,
            const std::shared_ptr<tom::ObservationAncilliarySimulationSettings>&
                ancilliarySettings,
            const tom::ObservableType observableType,
            const tom::LinkDefinition& linkEnds) {
            tss::addAncilliarySettingsToObservationSimulationSettings<
                TIME_TYPE, const tom::ObservableType,
                const tom::LinkDefinition&>(observationSimulationSettings,
                                            ancilliarySettings, observableType,
                                            linkEnds);
        }


    }  // namespace simulation_setup

}  // namespace tudat

namespace tudatpy {
    namespace numerical_simulation {
        namespace estimation_setup {
            namespace observation {

                PYBIND11_MODULE(expose_observation, m) {
                    py::module_::import(
                        "tudatpy.numerical_simulation.estimation");
                    py::module_::import("tudatpy.data");
                    // py::module_::import("tudatpy.data.io.expose_io");
                    // ################      Link Definition ################

                    py::enum_<tom::LinkEndType>(
                        m, "LinkEndType",
                        R"doc(Enumeration of available link end types.


	:member unidentified_link_end:
	:member transmitter:
	:member reflector1:
	:member retransmitter:
	:member reflector2:
	:member reflector3:
	:member reflector4:
	:member receiver:
	:member observed_body:
)doc")
                        .value("unidentified_link_end",
                               tom::LinkEndType::unidentified_link_end)
                        .value("transmitter", tom::LinkEndType::transmitter)
                        .value("reflector1", tom::LinkEndType::reflector1)
                        .value("retransmitter", tom::LinkEndType::retransmitter)
                        .value("reflector2", tom::LinkEndType::reflector2)
                        .value("reflector3", tom::LinkEndType::reflector3)
                        .value("reflector4", tom::LinkEndType::reflector4)
                        .value("receiver", tom::LinkEndType::receiver)
                        .value("observed_body", tom::LinkEndType::observed_body)
                        .export_values();


                    m.def(
                        "one_way_downlink_link_ends",
                        &tom::getOneWayDownlinkLinkEndsList,
                        py::arg("transmitter"), py::arg("receivers"),
                        R"doc(Function for defining one-way downlinks via LinkDefinition types.

	Function for defining single or multiple one-way downlinks.
	Multiple downlinks share the same transmitters, but may each have a different receiver.
	For each downlink, the returned list will contain an additional `LinkDefinition` type.


	:param transmitter:
		`LinkEndId` type (tuple of strings), where the first entrance identifies the body and the second entry the reference point of the single transmitter link end.

	:param receivers:
		List of `LinkEndId` types (tuple of strings), where for each tuple the first entrance identifies the body and the second entry the reference point of the receiver link end(s).

	:return:
		List of one or more `LinkDefinition` types, each defining the geometry for one one-way downlink.
		A `LinkDefinition` type for a one one-way link is composed a dict with one `receiver` and one `transmitter` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` key, to each of which a `LinkEndId` type is mapped.

)doc");

                    m.def(
                        "one_way_uplink_link_ends",
                        &tom::getOneWayUplinkLinkEndsList,
                        py::arg("transmitters"), py::arg("receiver"),
                        R"doc(Function for defining one-way uplinks via LinkDefinition types.

	Function for defining single or multiple one-way uplinks.
	Multiple uplinks share the same receiver, but may each have a different transmitter.
	For each uplink, the returned list will contain an additional `LinkDefinition` type.


	:param transmitters:
		List of `LinkEndId` types (tuple of strings), where for each tuple the first entrance identifies the body and the second entry the reference point of the transmitter link end(s).

	:param receiver:
		`LinkEndId` type (tuple of strings), where the first entrance identifies the body and the second entry the reference point of the single receiver link end.

	:return:
		List of one or more `LinkDefinition` types, each defining the geometry for one one-way uplink.
		A `LinkDefinition` type for a one one-way link is composed a dict with one `receiver` and one `transmitter` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` key, to each of which a `LinkEndId` type is mapped.

)doc");

                    m.def("get_default_reference_link_end",
                          &tom::getDefaultReferenceLinkEndType,
                          py::arg("observabl_type"), "");

                    // ###########      Observation Model Settings
                    // ################


                    py::class_<tom::LinkEndId, std::shared_ptr<tom::LinkEndId>>(
                        m, "LinkEndId",
                        R"doc(Object serving as identifier of a specific link end.

)doc")
                        .def_property_readonly("body_name",
                                               &tom::LinkEndId::getBodyName, "")
                        .def_property_readonly("reference_point",
                                               &tom::LinkEndId::getStationName,
                                               "");


                    m.def(
                        "body_origin_link_end_id",
                        py::overload_cast<const std::string&>(&tom::linkEndId),
                        py::arg("body_name"),
                        R"doc(Function to create a link end identifier for the origin (typically center of mass) of a body.

	Function to create a link end identifier for the origin (typically center of mass) of a body.
	Using this option will simulate the origin of a body transmitter, receiving, etc. the observation.
	Although this is typically not physically realistic, it can be a useful approximation, in particular for simulation studies.


	:param body_name:
		Name of the body

	:return:
		A LinkEndId object representing the center of mass of a body

)doc");

                    m.def(
                        "body_reference_point_link_end_id",
                        py::overload_cast<const std::string&,
                                          const std::string&>(&tom::linkEndId),
                        py::arg("body_name"), py::arg("reference_point_id"),
                        R"doc(Function to create a link end identifier for a reference point on a body.

	Function to create a link end identifier for a reference point on a body, where the reference point
	is typically the identifier of a ground stations


	:param body_name:
		Name of the body on which the reference point is located

	:param body_name:
		Name of the reference point on the body.

	:return:
		A LinkEndId object representing a reference point on a body

)doc");


                    py::class_<tom::LinkDefinition,
                               std::shared_ptr<tom::LinkDefinition>>(
                        m, "LinkDefinition",
                        R"doc(Object storing the link ends involved in a given observation.

)doc")
                        .def(py::init<const std::map<tom::LinkEndType,
                                                     tom::LinkEndId>&>(),
                             py::arg("link_ends"));
                    //            .def_property( "link_ends",
                    //            &tom::LinkDefinition::linkEnds_,
                    //                           );


                    m.def("link_definition", &tom::linkDefinition,
                          py::arg("link_ends"),
                          R"doc(Function to create a link definition object.

	:param link_ends:
		Dictionary of link ends, with the key denoting the role in the observaton, and the associated value the identifier for the link end.
	:return:
		The ``LinkDefinition`` object storing the link ends of the observation
)doc");


                    py::enum_<tom::ObservableType>(
                        m, "ObservableType",
                        R"doc(Enumeration of available observable types.


	:member one_way_range_type:
	:member n_way_range_type:
	:member angular_position_type:
	:member relative_angular_position_type:
	:member position_observable_type:
	:member velocity_observable_type:
	:member one_way_instantaneous_doppler_type:
	:member one_way_averaged_doppler_type:
	:member two_way_instantaneous_doppler_type:
	:member n_way_averaged_doppler_type:
	:member euler_angle_313_observable_type:
)doc")
                        .value("one_way_range_type",
                               tom::ObservableType::one_way_range)
                        .value("n_way_range_type",
                               tom::ObservableType::n_way_range)
                        .value("angular_position_type",
                               tom::ObservableType::angular_position)
                        .value("relative_angular_position_type",
                               tom::ObservableType::angular_position)
                        .value("position_observable_type",
                               tom::ObservableType::position_observable)
                        .value("velocity_observable_type",
                               tom::ObservableType::velocity_observable)
                        .value("one_way_instantaneous_doppler_type",
                               tom::ObservableType::one_way_doppler)
                        .value("one_way_averaged_doppler_type",
                               tom::ObservableType::one_way_differenced_range)
                        .value("two_way_instantaneous_doppler_type",
                               tom::ObservableType::two_way_doppler)
                        .value("n_way_averaged_doppler_type",
                               tom::ObservableType::n_way_differenced_range)
                        .value("euler_angle_313_observable_type",
                               tom::ObservableType::euler_angle_313_observable)
                        .value(
                            "dsn_one_way_averaged_doppler",
                            tom::ObservableType::dsn_one_way_averaged_doppler)
                        .value("dsn_n_way_averaged_doppler",
                               tom::ObservableType::dsn_n_way_averaged_doppler)
                        .export_values();


                    py::class_<
                        tom::DopplerProperTimeRateSettings,
                        std::shared_ptr<tom::DopplerProperTimeRateSettings>>(
                        m, "DopplerProperTimeRateSettings",
                        R"doc(Base class to defining proper time rate settings.

	Functional (base) class for settings of proper time rate (at a single link end) for instantaneous Doppler observation model settings.
	Specific proper time rate settings must be defined using an object derived from this class.
	The derived classes are made accessible via dedicated factory functions.

)doc");

                    py::class_<tom::ObservationModelSettings,
                               std::shared_ptr<tom::ObservationModelSettings>>(
                        m, "ObservationSettings",
                        R"doc(Base class for defining observation settings.

	Functional (base) class for settings of observation models.
	Observation model settings define at least the type and geometry of a given observation.
	They can furthermore set observation biases and/or light-time corrections.
	Simple observation models settings that are fully characterised by these elements can be managed by this base class, which can be instantiated through dedicated factory functions, such as
	:func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_range`, :func:`~tudatpy.numerical_simulation.estimation_setup.observation.cartesian_position`,
	:func:`~tudatpy.numerical_simulation.estimation_setup.observation.angular_position`, etc.
	Model settings for specific observation models that require additional information such as integration time, retransmission time, etc. must be defined using an object derived from this class.
	The derived classes are made accessible through further factory functions.

)doc");

                    py::class_<
                        tom::OneWayDopplerObservationSettings,
                        std::shared_ptr<tom::OneWayDopplerObservationSettings>,
                        tom::ObservationModelSettings>(
                        m, "OneWayDopplerObservationSettings",
                        R"doc(Class for defining the settings of one-way instantaneous Doppler observation models.

	:class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` derived class for one-way instantaneous Doppler observation model settings.
	Settings object can account for additional observation model aspects such as light time corrections and proper time rate settings.
	Instances of this class can be created via the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_doppler_instantaneous` factory function.

)doc");

                    py::class_<
                        tom::NWayRangeObservationSettings,
                        std::shared_ptr<tom::NWayRangeObservationSettings>,
                        tom::ObservationModelSettings>(
                        m, "NWayRangeObservationSettings", "");

                    py::class_<
                        tom::LightTimeConvergenceCriteria,
                        std::shared_ptr<tom::LightTimeConvergenceCriteria>>(
                        m, "LightTimeConvergenceCriteria", "");

                    py::enum_<tom::LightTimeFailureHandling>(
                        m, "LightTimeFailureHandling",
                        R"doc(Enumeration of behaviour when failing to converge light-time with required settings.


	:member accept_without_warning:
	:member print_warning_and_accept:
	:member throw_exception:
)doc")
                        .value("accept_without_warning",
                               tom::LightTimeFailureHandling::
                                   accept_without_warning)
                        .value("print_warning_and_accept",
                               tom::LightTimeFailureHandling::
                                   print_warning_and_accept)
                        .value("throw_exception",
                               tom::LightTimeFailureHandling::throw_exception)
                        .export_values();


                    m.def(
                        "light_time_convergence_settings",
                        &tom::lightTimeConvergenceCriteria,
                        py::arg("iterate_corrections") = false,
                        py::arg("maximum_number_of_iterations") = 50,
                        py::arg("absolute_tolerance") = TUDAT_NAN,
                        py::arg("failure_handling") =
                            tom::accept_without_warning,
                        R"doc(Factory function for creating settings for a one-way range observable.

	Factory function for creating observation model settings of one-way range type observables, for a single link definition. The associated observation model creates
	a single-valued observable :math:`h_{_{\text{1-range}}}` as follows (in the unbiased case):

	.. math::
	   h_{_{\text{1-range}}}(t_{R},t_{T})=|\mathbf{r}_{R}(t_{R})-\mathbf{r}_{T}(t_{T})| + \Delta s

	where :math:`\mathbf{r}_{R}`, :math:`\mathbf{r}_{T}`, :math:`t_{R}` and :math:`t_{T}` denote the position function of receiver and transmitter, and evaluation time
	of receiver and transmitter. The term :math:`\Delta s` denotes light-time corrections due to e.g relativistic, atmospheric effects (as defined by the ``light_time_correction_settings`` input).
	The transmission and reception time are related to the light-time :math:`T=t_{R}-t_{T}`, which is in turn related to the one-way range as :math:`T=h/c`
	As a result, the calculation of the one-way range (and light-time) requires the iterative solution of the equation:

	.. math::
	   t_{R}-t_{T}=c\left(|\mathbf{r}_{R}(t_{R})-\mathbf{r}(t_{R})| + \Delta s\right)

	 The method for the iterative solution is described in the :func:`light_time_convergence_settings` entry


	:param link_ends:
		Set of link ends that define the geometry of the observation. This observable requires the
		`transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.

	:param light_time_correction_settings:
		List of corrections for the light-time that are to be used. Default is none, which will result
		in the signal being modelled as moving in a straight line with the speed of light

	:param bias_settings:
		Settings for the observation bias that is to be used for the observation, default is None (unbiased observation)

	:param light_time_convergence_settings:
		Settings for convergence of the light-time (default settings defined in :func:`light_time_convergence_settings`)

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` class defining the settings for the one-way observable.
)doc");

                    m.def(
                        "one_way_range", &tom::oneWayRangeSettings,
                        py::arg("link_ends"),
                        py::arg("light_time_correction_settings") =
                            std::vector<std::shared_ptr<
                                tom::LightTimeCorrectionSettings>>(),
                        py::arg("bias_settings") = nullptr,
                        py::arg("light_time_convergence_settings") =
                            std::make_shared<
                                tom::LightTimeConvergenceCriteria>(),
                        R"doc(Factory function for creating settings for a one-way range observable.

	Factory function for creating observation model settings of one-way range type observables, for a single link definition. The associated observation model creates
	a single-valued observable :math:`h_{_{\text{1-range}}}` as follows (in the unbiased case):

	.. math::
	   h_{_{\text{1-range}}}(t_{R},t_{T})=|\mathbf{r}_{R}(t_{R})-\mathbf{r}_{T}(t_{T})| + \Delta s

	where :math:`\mathbf{r}_{R}`, :math:`\mathbf{r}_{T}`, :math:`t_{R}` and :math:`t_{T}` denote the position function of receiver and transmitter, and evaluation time
	of receiver and transmitter. The term :math:`\Delta s` denotes light-time corrections due to e.g relativistic, atmospheric effects (as defined by the ``light_time_correction_settings`` input).
	The transmission and reception time are related to the light-time :math:`T=t_{R}-t_{T}`, which is in turn related to the one-way range as :math:`T=h/c`
	As a result, the calculation of the one-way range (and light-time) requires the iterative solution of the equation:

	.. math::
	   t_{R}-t_{T}=c\left(|\mathbf{r}_{R}(t_{R})-\mathbf{r}(t_{R})| + \Delta s\right)

	 The method for the iterative solution is described in the :func:`light_time_convergence_settings` entry


	:param link_ends:
		Set of link ends that define the geometry of the observation. This observable requires the
		`transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.

	:param light_time_correction_settings:
		List of corrections for the light-time that are to be used. Default is none, which will result
		in the signal being modelled as moving in a straight line with the speed of light

	:param bias_settings:
		Settings for the observation bias that is to be used for the observation, default is None (unbiased observation)

	:param light_time_convergence_settings:
		Settings for convergence of the light-time (default settings defined in :func:`light_time_convergence_settings`)

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` class defining the settings for the one-way observable.
)doc");

                    m.def(
                        "two_way_range", &tom::twoWayRangeSimple,
                        py::arg("link_ends"),
                        py::arg("light_time_correction_settings") =
                            std::vector<std::shared_ptr<
                                tom::LightTimeCorrectionSettings>>(),
                        py::arg("bias_settings") = nullptr,
                        py::arg("light_time_convergence_settings") =
                            std::make_shared<
                                tom::LightTimeConvergenceCriteria>(),
                        R"doc(Factory function for creating settings for a two-way range observable.

	Same as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_range`, with :math:`n=2`. This function is provided
	for convenience.


	:param link_ends:
		Set of link ends that define the geometry of the observation. This observable requires the
		`transmitter`, `retransmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined

	:param light_time_correction_settings:
		List of corrections for the light-time that are to be used for each constituent one-way range. Default is none, which will result
		in the signal being modelled as moving in a straight line with the speed of light

	:param bias_settings:
		Settings for the observation bias that is to be used for the observation, default is none (unbiased observation).
		Note that only one bias setting is applied to the n-way observable.

	:param light_time_convergence_settings:
		Settings for convergence of the light-time (default settings defined in :func:`light_time_convergence_settings`)

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived :class:`NWayRangeObservationSettings` class.
)doc");

                    m.def(
                        "two_way_range_from_one_way_links", &tom::twoWayRange,
                        py::arg("one_way_range_settings"),
                        py::arg("bias_settings") = nullptr,
                        R"doc(Factory function for creating settings for a two-way range observable.

	Same as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_range_from_one_way_links`, with :math:`n=2`. This function is provided
	for convenience.


	:param one_way_range_settings:
		List of observation model settings of size two, with the first entry the one-way range settings for the uplink, and the second entry the one-way range settings for the downlink.
		The ``LinkDefinition`` of this two-way range observable is created from this list, with the ``transmitter`` and ``retransmitter1`` defined by the
		``transmitter`` and ``receiver`` of the first entry in this list. The ``retransmitter`` and ``receiver`` are defined by the
		``transmitter`` and ``receiver`` of the second entry of this list.

	:param bias_settings:
		Settings for the observation bias that is to be used for the observation, default is none (unbiased observation).
		Note that only one bias setting is applied to the n-way observable.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived :class:`NWayRangeObservationSettings` class.
)doc");

                    m.def(
                        "n_way_range", &tom::nWayRangeSimple,
                        py::arg("link_ends"),
                        py::arg("light_time_correction_settings") =
                            std::vector<std::shared_ptr<
                                tom::LightTimeCorrectionSettings>>(),
                        py::arg("bias_settings") = nullptr,
                        py::arg("light_time_convergence_settings") =
                            std::make_shared<
                                tom::LightTimeConvergenceCriteria>(),
                        R"doc(Factory function for creating settings for a n-way range observable.

	Factory function for creating observation model settings of n-way range type observables, for a single link definition. The associated observation model creates
	a single-valued observable :math:`h_{_{\text{N-range}}}` by combining together a series :math:`n` one-way range observations
	(see :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_range`). By default, the reception time of the :math:`i^{th}` one-way range is set as the
	transmission time of the :math:`(i+1)^{th}` one-way range. A retransmission delay may be defined by ancilliary settings (see TODO) when creating observation
	simulation setings.

	For this factory function, the settings for each constituent one-way range (with the exception of the link end identifiers) are equal.


	:param link_ends:
		Set of link ends that define the geometry of the observation. This observable requires the
		`transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined, as well
		as a `retransmitter1`, `retransmitter2`, .... (with the number of retransmitters to be defined by the user).

	:param light_time_correction_settings:
		List of corrections for the light-time that are to be used for each constituent one-way range. Default is none, which will result
		in the signal being modelled as moving in a straight line with the speed of light

	:param bias_settings:
		Settings for the observation bias that is to be used for the observation, default is none (unbiased observation).
		Note that only one bias setting is applied to the n-way observable.

	:param light_time_convergence_settings:
		Settings for convergence of the light-time (default settings defined in :func:`light_time_convergence_settings`)

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived :class:`NWayRangeObservationSettings` class.
)doc");

                    m.def(
                        "n_way_range_from_one_way_links", &tom::nWayRange,
                        py::arg("one_way_range_settings"),
                        py::arg("bias_settings") = nullptr,
                        R"doc(Factory function for creating settings for a n-way range observable.

	Factory function for creating observation model settings of n-way range type observables, for a single link definition. The
	implementation is the same as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_range`, with the difference
	that the constituent one-way ranges may have different settings.s


	:param one_way_range_settings:
		List of observation model settings for each of the :math:`n` constituent one-way ranges of the n-way range observable.
		The ``LinkDefinition`` of this n-way range observable is created from this list, with the ``transmitter`` and ``retransmitter1`` defined by the
		``transmitter`` and ``receiver`` of the first entry in this list. The ``retransmitter``(n-1) and ``receiver`` are defined by the
		``transmitter`` and ``receiver`` of the :math:`n`^{th} entry of this list.

	:param bias_settings:
		Settings for the observation bias that is to be used for the observation, default is none (unbiased observation).
		Note that only one bias setting is applied to the n-way observable.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived :class:`NWayRangeObservationSettings` class.
)doc");

                    m.def(
                        "angular_position", &tom::angularPositionSettings,
                        py::arg("link_ends"),
                        py::arg("light_time_correction_settings") =
                            std::vector<std::shared_ptr<
                                tom::LightTimeCorrectionSettings>>(),
                        py::arg("bias_settings") = nullptr,
                        py::arg("light_time_convergence_settings") =
                            std::make_shared<
                                tom::LightTimeConvergenceCriteria>(),
                        R"doc(Factory function for creating settings for an angular position observable.

	Factory function for creating observation model settings of angular position type observables (as right ascension :math:`\alpha` and declination :math:`\delta`),
	for a single link definition. The associated observation model creates an observable :math:`\mathbf{h}_{_{\text{ang.pos.}}}` of type two as follows (in the unbiased case):

	.. math::
	   \Delta\mathbf{r}=\mathbf{r}_{R}(t_{R})-\mathbf{r}_{T}(t_{T})\\
	   \tan\alpha=\frac{\Delta r_{y}}{\Delta r_{x}}\\
	   \delta=\frac{\Delta r_{z}}{\sqrt{\Delta r_{x}^{2}+\Delta r_{y}^{2}}}\\
	   \mathbf{h}_{_{\text{ang.pos.}}} = [\alpha;\delta]

	The relative position vector :math:`\Delta\mathbf{r}` is computed identically as described for the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_range`
	The angular position observable can be used for optical astrometry, VLBI, etc. Due to the definition of this observable, the xy-plane is defined by the global frame orientation of the
	environment.


	:param link_ends:
		Set of link ends that define the geometry of the observation. This observable requires the
		`transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.

	:param light_time_correction_settings:
		List of corrections for the light-time that are to be used. Default is none, which will result
		in the signal being modelled as moving in a straight line with the speed of light

	:param bias_settings:
		Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

	:param light_time_convergence_settings:
		Settings for convergence of the light-time (default settings defined in :func:`light_time_convergence_settings`)

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` class defining the settings for the angular position observable.
)doc");

                    m.def(
                        "relative_angular_position",
                        &tom::relativeAngularPositionSettings,
                        py::arg("link_ends"),
                        py::arg("light_time_correction_settings") =
                            std::vector<std::shared_ptr<
                                tom::LightTimeCorrectionSettings>>(),
                        py::arg("bias_settings") = nullptr,
                        py::arg("light_time_convergence_settings") =
                            std::make_shared<
                                tom::LightTimeConvergenceCriteria>(),
                        R"doc(Factory function for creating settings for an angular position observable.

	Factory function for creating observation model settings of angular position type observables (as right ascension :math:`\alpha` and declination :math:`\delta`),
	for a single link definition. The associated observation model creates an observable :math:`\mathbf{h}_{_{\text{ang.pos.}}}` of type two as follows (in the unbiased case):

	.. math::
	   \Delta\mathbf{r}=\mathbf{r}_{R}(t_{R})-\mathbf{r}_{T}(t_{T})\\
	   \tan\alpha=\frac{\Delta r_{y}}{\Delta r_{x}}\\
	   \delta=\frac{\Delta r_{z}}{\sqrt{\Delta r_{x}^{2}+\Delta r_{y}^{2}}}\\
	   \mathbf{h}_{_{\text{ang.pos.}}} = [\alpha;\delta]

	The relative position vector :math:`\Delta\mathbf{r}` is computed identically as described for the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_range`
	The angular position observable can be used for optical astrometry, VLBI, etc. Due to the definition of this observable, the xy-plane is defined by the global frame orientation of the
	environment.


	:param link_ends:
		Set of link ends that define the geometry of the observation. This observable requires the
		`transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.

	:param light_time_correction_settings:
		List of corrections for the light-time that are to be used. Default is none, which will result
		in the signal being modelled as moving in a straight line with the speed of light

	:param bias_settings:
		Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

	:param light_time_convergence_settings:
		Settings for convergence of the light-time (default settings defined in :func:`light_time_convergence_settings`)

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` class defining the settings for the angular position observable.
)doc");

                    m.def(
                        "cartesian_position", &tom::positionObservableSettings,
                        py::arg("link_ends"),
                        py::arg("bias_settings") = nullptr,
                        R"doc(Factory function for creating settings for a Cartesian position observable.

	Factory function for creating observation model settings of Cartesian position type observables.
	Note that this observable is typically not realized in reality, but can be very useful for verification or analysis purposes.
	This observable provides the inertial (w.r.t. global frame origin) Cartesian position of the `observed_body` defined by the `link_ends` input.
	The observable has size 3, and contains the :math:`x`, :math:`y` and :math:`z` position


	:param link_ends:
		Set of link ends that define the geometry of the observation. This observable requires that the
		`observed_body`` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.

	:param bias_settings:
		Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` class defining the settings for the cartesian position observable.
)doc");

                    m.def(
                        "relative_cartesian_position",
                        &tom::relativePositionObservableSettings,
                        py::arg("link_ends"),
                        py::arg("bias_settings") = nullptr,
                        R"doc(Factory function for creating settings for a Cartesian position observable.

	Factory function for creating observation model settings of Cartesian position type observables.
	Note that this observable is typically not realized in reality, but can be very useful for verification or analysis purposes.
	This observable provides the inertial (w.r.t. global frame origin) Cartesian position of the `observed_body` defined by the `link_ends` input.
	The observable has size 3, and contains the :math:`x`, :math:`y` and :math:`z` position


	:param link_ends:
		Set of link ends that define the geometry of the observation. This observable requires that the
		`observed_body`` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.

	:param bias_settings:
		Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` class defining the settings for the cartesian position observable.
)doc");


                    m.def(
                        "cartesian_velocity", &tom::velocityObservableSettings,
                        py::arg("link_ends"),
                        py::arg("bias_settings") = nullptr,
                        R"doc(Factory function for creating settings for a Cartesian velocity observable.

	Factory function for creating observation model settings of Cartesian position type observables.
	Note that this observable is typically not realized in reality, but can be very useful for verification or analysis purposes.
	This observable provides the inertial (w.r.t. global frame origin) Cartesian velocity of the `observed_body` defined by the `link_ends` input.
	The observable has size 3, and contains the :math:`x`, :math:`y` and :math:`z` velocity


	:param link_ends:
		Set of link ends that define the geometry of the observation. This observable requires that the
		`observed_body`` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.

	:param bias_settings:
		Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` class defining the settings for the cartesian velocity observable.
)doc");

                    m.def("euler_angles_313",
                          &tom::eulerAngle313ObservableSettings,
                          py::arg("link_ends"),
                          py::arg("bias_settings") = nullptr, "");


                    m.def(
                        "one_way_doppler_instantaneous",
                        &tom::oneWayOpenLoopDoppler, py::arg("link_ends"),
                        py::arg("light_time_correction_settings") =
                            std::vector<std::shared_ptr<
                                tom::LightTimeCorrectionSettings>>(),
                        py::arg("bias_settings") = nullptr,
                        py::arg("transmitter_proper_time_rate_settings") =
                            nullptr,
                        py::arg("receiver_proper_time_rate_settings") = nullptr,
                        py::arg("light_time_convergence_settings") =
                            std::make_shared<
                                tom::LightTimeConvergenceCriteria>(),
                        py::arg("normalized_with_speed_of_light") = false,
                        R"doc(Factory function for creating settings for a one-way instantaneous Doppler observable.

	Factory function for creating settings for a one-way instantaneous Doppler observable for a single link definition. The associated observation model creates
	a single-valued observable :math:`h_{_{\text{1-Dopp.}}}` as follows (in the unbiased case):

	.. math::
	   h_{_{\text{1-Dopp.}}}=c\left(\frac{d\tau_{T}}{dt_{T}}\frac{t_{T}}{dt_{R}}\frac{dt_{R}}{d\tau_{R}}-1\right)

	where :math:`t` and :math:`\tau` denote coordinate and proper time of the transmitter T and receiver R, respectively.
	The receiver and transmitter position and coordinate time are computed identically as described for the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_range`.
	The detailed mathematical implementation are described on TODO.

	This observable represents the 'instantaneous (non-integrated)' Doppler observable, as obtained from open-loop observations.
	It should *not* be used for the modelling of the typical closed-loop observations used in deep space tracking (for which the
	:func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_doppler_averaged` should be used)

	The coordinate
	time derivative :math:`\frac{t_{A}}{dt_{B}}` is always computed when generating this observable. Settings for the proper time
	rates :math:`\frac{d\tau}{dt}` can be specified by the user through the ``transmitter_proper_time_rate_settings`` and ``receiver_proper_time_rate_settings``
	inputs. Whene these are left empty, the proper time rates are omitted (set to 1.0).

	The observable may be non-dimensionalized by the speed of light :math:`c`, which results in the observable being equal to thee received and transmitted signal frequencies :math:`f_{R}/f_{T}-1`.


	:param link_ends:
		Set of link ends that define the geometry of the observation. This observable requires that the
		`transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.

	:param light_time_correction_settings:
		List of corrections for the light-time that are to be used. Default is none, which will result
		in the signal being modelled as moving in a straight line with the speed of light

	:param bias_settings:
		Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

	:param transmitter_proper_time_rate_settings:
		Settings for computing the transmitter proper time rate :math:`\frac{d\tau}{dt}`, default is none (:math:`\frac{d\tau}{dt}=1`)

	:param receiver_proper_time_rate_settings:
		Settings for computing the receiver proper time rate :math:`\frac{d\tau}{dt}`, default is none (:math:`\frac{d\tau}{dt}=1`)

	:param light_time_convergence_settings:
		Settings for convergence of the light-time (default settings defined in :func:`light_time_convergence_settings`)

	:param normalized_with_speed_of_light:
		Option to non-dimensionalize the observable with speed of light :math:`c`

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived :class:`OneWayDopplerObservationSettings` class defining the settings for the one-way open doppler observable observable.
)doc");

                    m.def(
                        "two_way_doppler_instantaneous_from_one_way_links",
                        py::overload_cast<
                            const std::shared_ptr<
                                tom::OneWayDopplerObservationSettings>,
                            const std::shared_ptr<
                                tom::OneWayDopplerObservationSettings>,
                            const std::shared_ptr<
                                tom::ObservationBiasSettings>>(
                            &tom::twoWayOpenLoopDoppler),
                        py::arg("uplink_doppler_settings"),
                        py::arg("downlink_doppler_settings"),
                        py::arg("bias_settings") = nullptr,
                        R"doc(Factory function for creating settings for a two-way instantaneous Doppler observable.


	Factory function for creating settings for a two-way instantaneous Doppler observable for a single link definition. The
	implementation is the same as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.two_way_doppler_instantaneous`, with the difference
	that the constituent one-way ranges may have different settings.

	The observable may be non-dimensionalized by the speed of light :math:`c` (in the constituent one-way Doppler observable settings),
	which results in the observable being equal to the received and transmitted signal frequencies :math:`f_{R}/f_{T}-1`.


	:param uplink_doppler_settings:
		Settings for uplink leg of one-way observable, created using :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_open_loop_doppler`

	:param downlink_doppler_settings:
		Settings for downlink leg of one-way observable, created using :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_open_loop_doppler`

	:param bias_settings:
		Settings for the observation bias that is to be used for the full observation, default is none (unbiased observation). Note that,
		even if no bias is applied to the two-way observable, the constituent one-way observables may still be biased.

	:param light_time_convergence_settings:
		Settings for convergence of the light-time (default settings defined in :func:`light_time_convergence_settings`)

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived :class:`TwoWayDopplerObservationSettings` class defining the settings for the two-way open doppler observable.
)doc");

                    m.def(
                        "two_doppler_instantaneous",
                        py::overload_cast<
                            const tom::LinkDefinition&,
                            const std::vector<std::shared_ptr<
                                tom::LightTimeCorrectionSettings>>&,
                            const std::shared_ptr<tom::ObservationBiasSettings>,
                            const std::shared_ptr<
                                tom::LightTimeConvergenceCriteria>,
                            const bool>(&tom::twoWayOpenLoopDoppler),
                        py::arg("link_ends"),
                        py::arg("light_time_correction_settings") =
                            std::vector<std::shared_ptr<
                                tom::LightTimeCorrectionSettings>>(),
                        py::arg("bias_settings") = nullptr,
                        py::arg("light_time_convergence_settings") =
                            std::make_shared<
                                tom::LightTimeConvergenceCriteria>(),
                        py::arg("normalized_with_speed_of_light") = false, "");

                    m.def(
                        "one_way_doppler_averaged",
                        py::overload_cast<
                            const tom::LinkDefinition&,
                            const std::vector<std::shared_ptr<
                                tom::LightTimeCorrectionSettings>>,
                            const std::shared_ptr<tom::ObservationBiasSettings>,
                            const std::shared_ptr<
                                tom::LightTimeConvergenceCriteria>>(
                            &tom::oneWayClosedLoopDoppler),
                        py::arg("link_ends"),
                        py::arg("light_time_correction_settings") =
                            std::vector<std::shared_ptr<
                                tom::LightTimeCorrectionSettings>>(),
                        py::arg("bias_settings") = nullptr,
                        py::arg("light_time_convergence_settings") =
                            std::make_shared<
                                tom::LightTimeConvergenceCriteria>(),
                        R"doc(Factory function for creating settings for a one-way averaged Doppler observable.

	Factory function for creating observation model settings for one-way averaged Doppler observables, for a single link definition. The associated observation model creates
	a single-valued observable :math:`h_{_{\text{1-\bar{Dopp}}}}` as follows (in the unbiased case):
	.. math::
	   h_{_{\text{1-\bar{Dopp}}}}&=c\int_{t-\Delta t}^{t+\Delta t}\frac{t_{T}}{dt_{R}}d\bar{t}\\
	                             &=\frac{h_{_{\text{1-range}}}(t_{R}=t+\Delta t,t_{T})-h_{_{\text{1-range}}}(t_{R}=t,t_{T})}{\Delta t}

	where, in the latter formulation (which is the one that is implemented), the observable is referenced to the receiver time. This averaged Doppler observable
	is computed as the difference of two one-way range observables (see :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_range`),
	with the reference time shifted by :math:`\Delta t`. As such, it is sensitive to numerical errors for small :math:`\Delta t`

	The integration time :math:`\Delta t` is defined in the ancilliary settings when simulating the observations (with 60 s as default).


	:param link_ends:
		Set of link ends that define the geometry of the observation. This observable requires that the
		`transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.

	:param light_time_correction_settings:
		List of corrections for the light-time that are to be used. Default is none, which will result
		in the signal being modelled as moving in a straight line with the speed of light

	:param bias_settings:
		Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

	:param light_time_convergence_settings:
		Settings for convergence of the light-time (default settings defined in :func:`light_time_convergence_settings`)

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived `OneWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable.
)doc");

                    m.def(
                        "two_way_doppler_averaged",
                        py::overload_cast<
                            const tom::LinkDefinition&,
                            const std::vector<std::shared_ptr<
                                tom::LightTimeCorrectionSettings>>,
                            const std::shared_ptr<tom::ObservationBiasSettings>,
                            const std::shared_ptr<
                                tom::LightTimeConvergenceCriteria>>(
                            &tom::twoWayDifferencedRangeObservationSettings),
                        py::arg("link_ends"),
                        py::arg("light_time_correction_settings") =
                            std::vector<std::shared_ptr<
                                tom::LightTimeCorrectionSettings>>(),
                        py::arg("bias_settings") = nullptr,
                        py::arg("light_time_convergence_settings") =
                            std::make_shared<
                                tom::LightTimeConvergenceCriteria>(),
                        R"doc(Factory function for creating settings for an n-way averaged Doppler observable.

	Factory function for creating observation model settings for n-way averaged Doppler observables, for a single link definition. The implemenation is
	analogous to the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_doppler_averaged` observable. But, in the present case
	the observable is computed from the difference of two n-way range observables, with the reference time shifted by :math:`\Delta t`.

	The integration time :math:`\Delta t` is defined in the ancilliary settings when simulating the observations (with 60 s as default).


	:param link_ends:
		Set of link ends that define the geometry of the observation. This observable requires the
		`transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined, as well
		as a `retransmitter1`, `retransmitter2`, .... (with the number of retransmitters to be defined by the user).

	:param light_time_correction_settings:
		List of corrections for the light-time that are to be used. Default is none, which will result
		in the signal being modelled as moving in a straight line with the speed of light

	:param bias_settings:
		Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

	:param light_time_convergence_settings:
		Settings for convergence of the light-time (default settings defined in :func:`light_time_convergence_settings`)

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived `NWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable.
)doc");

                    m.def(
                        "two_way_doppler_averaged_from_one_way_links",
                        py::overload_cast<const std::vector<std::shared_ptr<
                                              tom::ObservationModelSettings>>,
                                          const std::shared_ptr<
                                              tom::ObservationBiasSettings>>(
                            &tom::twoWayDifferencedRangeObservationSettings),
                        py::arg("one_way_range_settings"),
                        py::arg("bias_settings") = nullptr,
                        R"doc(Factory function for creating settings for an n-way averaged Doppler observable.

	Factory function for creating observation model settings for n-way averaged Doppler observables, for a single link definition. The implemenation is
	analogous to the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_doppler_averaged` observable. But, in the present case
	the observable is computed from the difference of two n-way range observables, with the reference time shifted by :math:`\Delta t`.

	The integration time :math:`\Delta t` is defined in the ancilliary settings when simulating the observations (with 60 s as default).


	:param link_ends:
		Set of link ends that define the geometry of the observation. This observable requires the
		`transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined, as well
		as a `retransmitter1`, `retransmitter2`, .... (with the number of retransmitters to be defined by the user).

	:param light_time_correction_settings:
		List of corrections for the light-time that are to be used. Default is none, which will result
		in the signal being modelled as moving in a straight line with the speed of light

	:param bias_settings:
		Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

	:param light_time_convergence_settings:
		Settings for convergence of the light-time (default settings defined in :func:`light_time_convergence_settings`)

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived `NWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable.
)doc");

                    m.def(
                        "n_way_doppler_averaged",
                        py::overload_cast<
                            const tom::LinkDefinition&,
                            const std::vector<std::shared_ptr<
                                tom::LightTimeCorrectionSettings>>,
                            const std::shared_ptr<tom::ObservationBiasSettings>,
                            const std::shared_ptr<
                                tom::LightTimeConvergenceCriteria>>(
                            &tom::nWayDifferencedRangeObservationSettings),
                        py::arg("link_ends"),
                        py::arg("light_time_correction_settings") =
                            std::vector<std::shared_ptr<
                                tom::LightTimeCorrectionSettings>>(),
                        py::arg("bias_settings") = nullptr,
                        py::arg("light_time_convergence_settings") =
                            std::make_shared<
                                tom::LightTimeConvergenceCriteria>(),
                        R"doc(Factory function for creating settings for an n-way averaged Doppler observable.

	Factory function for creating observation model settings for n-way averaged Doppler observables, for a single link definition. The implemenation is
	analogous to the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_doppler_averaged` observable. But, in the present case
	the observable is computed from the difference of two n-way range observables, with the reference time shifted by :math:`\Delta t`.

	The integration time :math:`\Delta t` is defined in the ancilliary settings when simulating the observations (with 60 s as default).


	:param link_ends:
		Set of link ends that define the geometry of the observation. This observable requires the
		`transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined, as well
		as a `retransmitter1`, `retransmitter2`, .... (with the number of retransmitters to be defined by the user).

	:param light_time_correction_settings:
		List of corrections for the light-time that are to be used. Default is none, which will result
		in the signal being modelled as moving in a straight line with the speed of light

	:param bias_settings:
		Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

	:param light_time_convergence_settings:
		Settings for convergence of the light-time (default settings defined in :func:`light_time_convergence_settings`)

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived `NWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable.
)doc");

                    m.def(
                        "n_way_doppler_averaged_from_one_way_links",
                        py::overload_cast<
                            const std::vector<
                                std::shared_ptr<tom::ObservationModelSettings>>,
                            const std::shared_ptr<tom::ObservationBiasSettings>,
                            const std::shared_ptr<
                                tom::LightTimeConvergenceCriteria>>(
                            &tom::nWayDifferencedRangeObservationSettings),
                        py::arg("one_way_range_settings"),
                        py::arg("bias_settings") = nullptr,
                        py::arg("light_time_convergence_settings") =
                            std::make_shared<
                                tom::LightTimeConvergenceCriteria>(),
                        R"doc(Factory function for creating settings for an n-way averaged Doppler observable.

	Factory function for creating observation model settings for n-way averaged Doppler observables, for a single link definition.
	The implementation is the same as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_doppler_averaged`, with the difference
	that the constituent one-way range observables may have different settings.


	:param one_way_range_settings:
		List of observation model settings for each of the :math:`n` constituent one-way ranges of the n-way averaged range rate observable.
		The ``LinkDefinition`` of this n-way range observable is created from this list, with the ``transmitter`` and ``retransmitter1`` defined by the
		``transmitter`` and ``receiver`` of the first entry in this list. The ``retransmitter``(n-1) and ``receiver`` are defined by the
		``transmitter`` and ``receiver`` of the :math:`n`^{th} entry of this list.

	:param bias_settings:
		Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived `NWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable.
)doc");

                    m.def(
                        "dsn_n_way_doppler_averaged",
                        py::overload_cast<
                            const tom::LinkDefinition&,
                            const std::vector<std::shared_ptr<
                                tom::LightTimeCorrectionSettings>>,
                            const std::shared_ptr<tom::ObservationBiasSettings>,
                            const std::shared_ptr<
                                tom::LightTimeConvergenceCriteria>>(
                            &tom::dsnNWayAveragedDopplerObservationSettings),
                        py::arg("link_ends"),
                        py::arg("light_time_correction_settings") =
                            std::vector<std::shared_ptr<
                                tom::LightTimeCorrectionSettings>>(),
                        py::arg("bias_settings") = nullptr,
                        py::arg("light_time_convergence_settings") =
                            std::make_shared<
                                tom::LightTimeConvergenceCriteria>(),
                        "");

                    m.def(
                        "dsn_n_way_doppler_averaged_from_one_way_links",
                        py::overload_cast<
                            const std::vector<
                                std::shared_ptr<tom::ObservationModelSettings>>,
                            const std::shared_ptr<tom::ObservationBiasSettings>,
                            const std::shared_ptr<
                                tom::LightTimeConvergenceCriteria>>(
                            &tom::dsnNWayAveragedDopplerObservationSettings),
                        py::arg("one_way_range_settings"),
                        py::arg("bias_settings") = nullptr,
                        py::arg("light_time_convergence_settings") =
                            std::make_shared<
                                tom::LightTimeConvergenceCriteria>(),
                        "");

                    py::class_<
                        tom::LightTimeCorrectionSettings,
                        std::shared_ptr<tom::LightTimeCorrectionSettings>>(
                        m, "LightTimeCorrectionSettings",
                        R"doc(Base class to defining light time correction settings.

	Functional (base) class for settings of light time corrections.
	This class is not used for calculations of corrections, but is used for the purpose of defining the light time correction properties.
	Specific light time correction settings must be defined using an object derived from this class.
	The derived classes are made accessible via dedicated factory functions, such as e.g.
	:func:`~tudatpy.numerical_simulation.estimation_setup.observation.first_order_relativistic_light_time_correction`

)doc");

                    m.def(
                        "first_order_relativistic_light_time_correction",
                        &tom::firstOrderRelativisticLightTimeCorrectionSettings,
                        py::arg("perturbing_bodies"),
                        R"doc(Factory function for creating settings for first-order relativistic light-time corrections.

	Factory function for creating settings for first-order relativistic light-time corrections: the correction to
	the light time of a (set of) stationary point masses, computed up to c−2 according to general relativity as formulated by e.g. Moyer (2000).
	One ambiguity in the model is the time at which the states of the perturbing bodies are evaluated. We distinguish two cases:

	* In the case where the perturbing body contains a link end of the observation (for instance perturbation due to Earth gravity field,
	  with one of the link ends being an Earth-based station), the time at which the Earth’s state is evaluated equals the transmission time if Earth acts as transmitter, and reception time if
	  Earth acts as receiver.
	* In other cases, where the perturbing body is not involved in the link ends, its state is evaluated at the midpoint time between transmitter and receiver.


	:param perturbing_bodies:
		A list containing the names of the bodies due to which the light-time correction is to be taken into account.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LightTimeCorrectionSettings` derived :class:`FirstOrderRelativisticLightTimeCorrectionSettings` class,
		defining the settings for the light-time corrections

)doc");

                    py::enum_<tom::TroposphericMappingModel>(
                        m, "TroposphericMappingModel", "")
                        .value("simplified_chao",
                               tom::TroposphericMappingModel::simplified_chao)
                        .value("niell", tom::TroposphericMappingModel::niell)
                        .export_values();

                    py::enum_<tom::WaterVaporPartialPressureModel>(
                        m, "WaterVaporPartialPressureModel", "")
                        .value("tabulated",
                               tom::WaterVaporPartialPressureModel::tabulated)
                        .value("bean_and_dutton",
                               tom::WaterVaporPartialPressureModel::
                                   bean_and_dutton)
                        .export_values();

                    m.def("dsn_tabulated_tropospheric_light_time_correction",
                          &tom::tabulatedTroposphericCorrectionSettings,
                          py::arg("file_names"),
                          py::arg("body_with_atmosphere_name") = "Earth",
                          py::arg("mapping_model") =
                              tom::TroposphericMappingModel::niell,
                          "");

                    m.def("saastamoinen_tropospheric_light_time_correction",
                          &tom::saastamoinenTroposphericCorrectionSettings,
                          py::arg("body_with_atmosphere_name") = "Earth",
                          py::arg("mapping_model") =
                              tom::TroposphericMappingModel::niell,
                          py::arg("water_vapor_partial_pressure_model") =
                              tom::WaterVaporPartialPressureModel::tabulated,
                          "");

                    m.def("dsn_tabulated_ionospheric_light_time_correction",
                          &tom::tabulatedIonosphericCorrectionSettings,
                          py::arg("file_names"),
                          py::arg("spacecraft_name_per_id") =
                              std::map<int, std::string>(),
                          py::arg("quasar_name_per_id") =
                              std::map<int, std::string>(),
                          py::arg("reference_frequency") = 2295e6,
                          py::arg("body_with_atmosphere_name") = "Earth", "");

                    m.def("jakowski_ionospheric_light_time_correction",
                          &tom::jakowskiIonosphericCorrectionSettings,
                          py::arg("ionosphere_height") = 400.0e3,
                          py::arg("first_order_delay_coefficient") = 40.3,
                          py::arg("solar_activity_data") = tudat::input_output::
                              solar_activity::readSolarActivityData(
                                  tudat::paths::getSpaceWeatherDataPath() +
                                  "/sw19571001.txt"),
                          py::arg("geomagnetic_pole_latitude") =
                              tuc::convertDegreesToRadians(80.9),
                          py::arg("geomagnetic_pole_longitude") =
                              tuc::convertDegreesToRadians(-72.6),
                          py::arg("use_utc_for_local_time_computation") = false,
                          py::arg("body_with_atmosphere_name") = "Earth", "");

                    m.def(
                        "inverse_power_series_solar_corona_light_time_"
                        "correction",
                        &tom::inversePowerSeriesSolarCoronaCorrectionSettings,
                        py::arg("coefficients") =
                            std::vector<double>{1.31 * 5.97e-6},
                        py::arg("positive_exponents") =
                            std::vector<double>{2.0},
                        py::arg("delay_coefficient") = 40.3,
                        py::arg("sun_body_name") = "Sun", "");

                    py::class_<tom::ObservationBiasSettings,
                               std::shared_ptr<tom::ObservationBiasSettings>>(
                        m, "ObservationBiasSettings",
                        R"doc(Base class to defining observation bias settings.

	Functional (base) class for settings of observation bias.
	Specific observation bias settings must be defined using an object derived from this class.
	The derived classes are made accessible via dedicated factory functions.

)doc");

                    m.def(
                        "absolute_bias", &tom::constantAbsoluteBias,
                        py::arg("bias_value"),
                        R"doc(Factory function for creating settings for an absolute observation bias.

	Factory function for creating settings for an absolute observation bias. When calculating the observable value, applying this setting
	will take the physically ideal observation :math:`h`, and modify it to obtain the biased observation :math:`\tilde{h}` as follows:

	.. math::
	   \tilde{h}=h+K

	where :math:`K` is the `bias_value`. For an observable with size greater than 1, :math:`K` is a vector and the addition is component-wise.


	:param bias_value:
		A vector containing the bias that is to be applied to the observable. This vector should be the same size as the observable to which it is
		applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`ConstantObservationBiasSettings` class, defining the settings for a constant, absolute observation bias.
)doc");

                    m.def(
                        "relative_bias", &tom::constantRelativeBias,
                        py::arg("bias_value"),
                        R"doc(Factory function for creating settings for a relative observation bias.

	Factory function for creating settings for a relative observation bias. When calculating the observable value, applying this setting
	will take the physically ideal observation :math:`h`, and modify it to obtain the biased observation :math:`\tilde{h}` as follows:

	.. math::
	   \tilde{h}=h(1+K)

	where :math:`K` is the`bias_value`. For an observable with size greater than 1, :math:`K` is a vector and the multiplication is component-wise.


	:param bias_value:
		A vector containing the bias that is to be applied to the observable. This vector should be the same size as the observable to which it is
		applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`ConstantObservationBiasSettings` class,
		defining the settings for a constant, relative observation bias.

)doc");

                    m.def(
                        "arcwise_absolute_bias",
                        py::overload_cast<const std::vector<double>&,
                                          const std::vector<Eigen::VectorXd>&,
                                          const tom::LinkEndType>(
                            &tom::arcWiseAbsoluteBias),
                        py::arg("arc_start_times"), py::arg("bias_values"),
                        py::arg("reference_link_end_type"),
                        R"doc(Factory function for creating settings for arc-wise absolute observation biases.

	Factory function for creating settings for arc-wise absolute observation biases.
	This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.absolute_bias` setting only through the option of setting the `bias_value` :math:`K` to a different values for each arc.


	:param bias_values_per_start_time:
		Dictionary, in which the bias value vectors for each arc are directly mapped to the starting times of the respective arc.
		The vectors should be the same size as the observable to which it is applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

	:param reference_link_end_type:
		Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.
)doc");

                    m.def(
                        "arcwise_absolute_bias",
                        py::overload_cast<
                            const std::map<double, Eigen::VectorXd>&,
                            const tom::LinkEndType>(&tom::arcWiseAbsoluteBias),
                        py::arg("bias_values_per_start_time"),
                        py::arg("reference_link_end_type"),
                        R"doc(Factory function for creating settings for arc-wise absolute observation biases.

	Factory function for creating settings for arc-wise absolute observation biases.
	This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.absolute_bias` setting only through the option of setting the `bias_value` :math:`K` to a different values for each arc.


	:param bias_values_per_start_time:
		Dictionary, in which the bias value vectors for each arc are directly mapped to the starting times of the respective arc.
		The vectors should be the same size as the observable to which it is applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

	:param reference_link_end_type:
		Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.
)doc");

                    m.def(
                        "arcwise_relative_bias",
                        py::overload_cast<const std::vector<double>&,
                                          const std::vector<Eigen::VectorXd>&,
                                          const tom::LinkEndType>(
                            &tom::arcWiseRelativeBias),
                        py::arg("arc_start_times"), py::arg("bias_values"),
                        py::arg("reference_link_end_type"),
                        R"doc(Factory function for creating settings for arc-wise relative observation biases.

	Factory function for creating settings for arc-wise relative observation biases.
	This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.relative_bias` setting only through the option of setting the `bias_value` :math:`K` to a different values for each arc.


	:param bias_values_per_start_time:
		Dictionary, in which the bias value vectors for each arc are directly mapped to the starting times of the respective arc.
		The vectors should be the same size as the observable to which it is applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

	:param reference_link_end_type:
		Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.
)doc");

                    m.def(
                        "arcwise_relative_bias",
                        py::overload_cast<
                            const std::map<double, Eigen::VectorXd>&,
                            const tom::LinkEndType>(&tom::arcWiseRelativeBias),
                        py::arg("bias_values_per_start_time"),
                        py::arg("reference_link_end_type"),
                        R"doc(Factory function for creating settings for arc-wise relative observation biases.

	Factory function for creating settings for arc-wise relative observation biases.
	This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.relative_bias` setting only through the option of setting the `bias_value` :math:`K` to a different values for each arc.


	:param bias_values_per_start_time:
		Dictionary, in which the bias value vectors for each arc are directly mapped to the starting times of the respective arc.
		The vectors should be the same size as the observable to which it is applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

	:param reference_link_end_type:
		Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.
)doc");

                    m.def("time_drift_bias", &tom::constantTimeDriftBias,
                          py::arg("bias_value"), py::arg("time_link_end"),
                          py::arg("ref_epoch"));

                    m.def("arc_wise_time_drift_bias",
                          py::overload_cast<const std::vector<Eigen::VectorXd>&,
                                            const std::vector<double>&,
                                            const tom::LinkEndType,
                                            const std::vector<double>&>(
                              &tom::arcWiseTimeDriftBias),
                          py::arg("bias_value"), py::arg("arc_start_times"),
                          py::arg("time_link_end"), py::arg("ref_epochs"));

                    m.def(
                        "arc_wise_time_drift_bias",
                        py::overload_cast<
                            const std::map<double, Eigen::VectorXd>&,
                            const tom::LinkEndType, const std::vector<double>>(
                            &tom::arcWiseTimeDriftBias),
                        py::arg("bias_value_per_start_time"),
                        py::arg("time_link_end"), py::arg("ref_epochs"));

                    m.def(
                        "combined_bias", &tom::multipleObservationBiasSettings,
                        py::arg("bias_list"),
                        R"doc(Factory function for creating settings for a combined observation bias.

	Factory function for creating settings for a combined observation bias, calculating by combining any number of bias types.
	Each contribution of the combined bias is computed from the unbiased observable, so when applying both a relative and absolute bias, we get:

	.. math::
	   \tilde{h}=h+K_{a}+hK_{r}

	And, crucially:

	.. math::
	   \tilde{h}\neq (h+K_{a})(1+K_{r})

	where :math:`K_{r}` and :math:`K_{a}` is the relative and absolute bias, respectively.


	:param bias_list:
		A list containing the bias the bias settings that are to be applied to the observable.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`MultipleObservationBiasSettings` class, combining the settings for multiple observation biases.
)doc");


                    // ###########    Observation Simulation Settings
                    // #############

                    py::enum_<tom::ObservationViabilityType>(
                        m, "ObservationViabilityType",
                        R"doc(Enumeration of observation viability criterion types.


	:member minimum_elevation_angle:
	:member body_avoidance_angle:
	:member body_occultation:
)doc")
                        .value("minimum_elevation_angle",
                               tom::ObservationViabilityType::
                                   minimum_elevation_angle)
                        .value(
                            "body_avoidance_angle",
                            tom::ObservationViabilityType::body_avoidance_angle)
                        .value("body_occultation",
                               tom::ObservationViabilityType::body_occultation)
                        .export_values();

                    py::enum_<tom::ObservationAncilliarySimulationVariable>(
                        m, "ObservationAncilliarySimulationVariable", "")
                        .value("link_ends_delays",
                               tom::ObservationAncilliarySimulationVariable::
                                   link_ends_delays)
                        .value("doppler_integration_time",
                               tom::ObservationAncilliarySimulationVariable::
                                   doppler_integration_time)
                        .value("doppler_reference_frequency",
                               tom::ObservationAncilliarySimulationVariable::
                                   doppler_reference_frequency)
                        .value("frequency_bands",
                               tom::ObservationAncilliarySimulationVariable::
                                   frequency_bands)
                        .value("reception_reference_frequency_band",
                               tom::ObservationAncilliarySimulationVariable::
                                   reception_reference_frequency_band)
                        .export_values();

                    py::class_<
                        tom::ObservationViabilitySettings,
                        std::shared_ptr<tom::ObservationViabilitySettings>>(
                        m, "ObservationViabilitySettings",
                        R"doc(Enumeration of observation viability criterion types.


	:member minimum_elevation_angle:
	:member body_avoidance_angle:
	:member body_occultation:
)doc");


                    m.def(
                        "elevation_angle_viability",
                        py::overload_cast<
                            const std::pair<std::string, std::string>,
                            const double>(
                            &tom::elevationAngleViabilitySettings),
                        py::arg("link_end_id"), py::arg("elevation_angle"),
                        R"doc(Factory function for defining single elevation angle viability setting.

	Factory function for defining elevation angle viability settings for single link end.
	When simulating observations, this setting ensures that any applicable observations, for which the local elevation angle at link end is less than some limit value, will be omitted.


	:param link_end_id:
		Link end (as defined by body/reference point pair, see TODO), for which the elevation angle viability setting is to be created.
		To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""].

	:param elevation_angle:
		Limit elevation angle, below which no observations are produced when using the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.simulate_observations` function. Note: this
		value must be in radians.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` class, defining the settings for observation viability
)doc");

                    m.def(
                        "body_avoidance_viability",
                        py::overload_cast<
                            const std::pair<std::string, std::string>,
                            const std::string, const double>(
                            &tom::bodyAvoidanceAngleViabilitySettings),
                        py::arg("link_end_id"), py::arg("body_to_avoid"),
                        py::arg("avoidance_angle"),
                        R"doc(Factory function for defining body avoidance observation viability settings.

	Factory function for defining body avoidance observation viability settings for single link ends.
	When simulating observations, this settings ensures that any applicable observations, for which the signal path passes 'too close' to a body, will be omitted.
	The definition of 'too close' is computed as the angle between:

	* The line-of-sight vector from a link end to a given third body
	* The line-of-sight between two link ends

	This constraint is typically used to prevent the Sun from being too close to the field-of-view of the telescope(s), as defined by
	a so-called 'SPE' (Sun-Probe-Earth) angle constraint. The present viability setting generalizes this constraint.


	:param link_end_id:
		Link end (as defined by body/reference point pair, see TODO), for which the viability settings are to be created.
		To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""] is entry in this list.
		For each link end included in this list, it will be checked if a signal received by and/or transmitted (or reflected) by this
		link end passes too close to the specified body.

	:param body_to_avoid:
		Name of the body which the signal path should not pass 'too close' to.

	:param avoidance_angle:
		Limit angle (generalization of SPE angle), below which no observations are produced when using the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.simulate_observations` function. Note: this
		value must be in radians.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings`, defining the settings for observation viability.
)doc");

                    m.def(
                        "body_occultation_viability",
                        py::overload_cast<
                            const std::pair<std::string, std::string>,
                            const std::string>(
                            &tom::bodyOccultationViabilitySettings),
                        py::arg("link_end_id"), py::arg("occulting_body"),
                        R"doc(Factory function for defining body occultation viability settings.

	Factory function for defining body occultation viability settings for single link ends.
	When simulating observations, this setting ensures that any applicable observations, for which the signal path is occulted by a given body, will be omitted.
	The occultation is computed using the shape model of the specified body.


	:param link_end_id:
		Link end (as defined by body/reference point pair, see TODO), for which the viability settings are to be created.
		To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""] is entry in this list.

	:param body_to_avoid:
		Name of the body which the signal path should not be occulted by.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings`, defining the settings for observation viability.
)doc");

                    m.def(
                        "elevation_angle_viability_list",
                        py::overload_cast<const std::vector<std::pair<
                                              std::string, std::string>>,
                                          const double>(
                            &tom::elevationAngleViabilitySettings),
                        py::arg("link_end_ids"), py::arg("elevation_angle"),
                        R"doc(Factory function for defining list of elevation angle viability settings.

	Factory function for defining elevation angle viability settings for multiple link ends.
	Each entry in the returned list contains the observation viability settings for one link end.
	When simulating observations, these settings ensure that any applicable observations, for which the local elevation angle at a link end is less than some limit value, will be omitted.


	:param link_end_ids:
		List of individual link ends (as defined by body/reference point pair, see TODO), for which the elevation angle viability setting is to be created.
		To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""].
		For each link end included in this list, it will be checked if a signal received by and/or transmitted (or reflected) by this
		link end violates the minimum elevation angle constraint.

	:param elevation_angle:
		Limit elevation angle, below which no observations are produced when using the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.simulate_observations` function. Note: this
		value must be in radians.

	:return:
		List of :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` objects, each defining the settings for observation viability of one link end.
)doc");

                    m.def(
                        "body_avoidance_viability_list",
                        py::overload_cast<const std::vector<std::pair<
                                              std::string, std::string>>,
                                          const std::string, const double>(
                            &tom::bodyAvoidanceAngleViabilitySettings),
                        py::arg("link_end_ids"), py::arg("body_to_avoid"),
                        py::arg("avoidance_angle"),
                        R"doc(Factory function for defining list of body avoidance viability settings.

	Factory function for defining body avoidance viability settings for multiple link ends.
	Each entry in the returned list contains the observation viability settings for one link end.
	When simulating observations, these settings ensure that any applicable observations, for which the signal path passes 'too close' to a body, will be omitted.
	The definition of 'too close' is computed as the angle between:

	* The line-of-sight vector from a link end to a given third body
	* The line-of-sight between two link ends

	This constraint is typically used to prevent the Sun from being too close to the field-of-view of the telescope(s), as defined by
	a so-called 'SPE' (Sun-Probe-Earth) angle constraint. The present viability setting generalizes this constraint.


	:param link_end_ids:
		List of individual link ends (as defined by body/reference point pair, see TODO), for which the elevation angle viability setting is to be created.
		To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""].

	:param body_to_avoid:
		Name of the body which the signal path should not pass 'too close' to.

	:param avoidance_angle:
		Limit angle (generalization of SPE angle), below which no observations are produced when using the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.simulate_observations` function. Note: this
		value must be in radians.

	:return:
		List of :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` objects, each defining the settings for observation viability of one link end.
)doc");

                    m.def(
                        "body_occultation_viability_list",
                        py::overload_cast<
                            const std::pair<std::string, std::string>,
                            const std::string>(
                            &tom::bodyOccultationViabilitySettings),
                        py::arg("link_end_id"), py::arg("occulting_body"),
                        R"doc(Factory function for defining body occultation viability settings.

	Factory function for defining body occultation viability settings for multiple link ends.
	Each entry in the returned list contains the observation viability settings for one link end.
	When simulating observations, these settings ensure that any applicable observations, for which the signal path is occulted by a given body, will be omitted.
	The occultation is computed using the shape model of the specified body.


	:param link_end_ids:
		List of individual link ends (as defined by body/reference point pair, see TODO), for which the viability settings are to be created.
		To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""] is entry in this list.
		For each link end included in this list, it will be checked if a signal received by and/or transmitted (or reflected) by this
		link end is occulted by the specified body.

	:param body_to_avoid:
		Name of the body which the signal path should not be occulted by.

	:return:
		List of :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` objects, each defining the settings for observation viability of one link end.
)doc");


                    py::class_<tss::ObservationSimulationSettings<double>,
                               std::shared_ptr<
                                   tss::ObservationSimulationSettings<double>>>(
                        m, "ObservationSimulationSettings",
                        R"doc(Base class for defining settings for simulating observations.

	Base class for defining settings for simulating observations.
	This simulation settings object defines observation times, noise and viability criteria, *etc.* at which observations are to be simulated.
	Therefore, one simulation settings object of this type can only refer to one combination of observable type and link geometry (LinkDefinition).
	The user does not interact with this class directly, but defines specific observation simulation settings using an object derived from this class (created through the associated factory function).

)doc")
                        .def_property("viability_settings_list",
                                      &tss::ObservationSimulationSettings<
                                          double>::getViabilitySettingsList,
                                      &tss::ObservationSimulationSettings<
                                          double>::setViabilitySettingsList,
                                      "")
                        .def_property(
                            "noise_function",
                            &tss::ObservationSimulationSettings<
                                double>::getObservationNoiseFunction,
                            py::overload_cast<
                                const std::function<double(const double)>&>(
                                &tss::ObservationSimulationSettings<
                                    double>::setObservationNoiseFunction),
                            "")
                        .def_property("observable_type",
                                      &tss::ObservationSimulationSettings<
                                          double>::getObservableType,
                                      &tss::ObservationSimulationSettings<
                                          double>::setObservableType,
                                      "")
                        .def_property_readonly(
                            "link_ends",
                            &tss::ObservationSimulationSettings<
                                double>::getLinkEnds,
                            "");


                    py::class_<
                        tss::TabulatedObservationSimulationSettings<double>,
                        std::shared_ptr<
                            tss::TabulatedObservationSimulationSettings<
                                double>>,
                        tss::ObservationSimulationSettings<double>>(
                        m, "TabulatedObservationSimulationSettings",
                        R"doc(Class for defining settings for simulating observations at a predefined set of times.

	:class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` derived class for defining settings for simulating observations
	at a predefined set of times
	This type defines predefined time epochs at which applicable observations are to be simulated, stored in a rigid, "tabulated" form.
	Some observation times may be discarded due to the use of viability settings.
	Instances of this class are typicall created via the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.tabulated_simulation_settings`
	and :func:`~tudatpy.numerical_simulation.estimation_setup.observation.tabulated_simulation_settings_list` factory functions.

)doc");


                    py::class_<
                        tom::ObservationAncilliarySimulationSettings,
                        std::shared_ptr<
                            tom::ObservationAncilliarySimulationSettings>>(
                        m, "ObservationAncilliarySimulationSettings",
                        R"doc(Class for holding ancilliary settings for observation simulation.

)doc")
                        .def("get_float_settings",
                             &tom::ObservationAncilliarySimulationSettings::
                                 getAncilliaryDoubleData,
                             py::arg("setting_type"),
                             py::arg("throw_exception") = true,
                             R"doc(
	:param setting_type:
		Type of the setting for which the value is to be returned

	:param throw_exception:
		Boolean defining whether to throw an exception if the requested setting does not exist, or does not exist as a floating point value.

	:return:
		Value of the requested ancilliary variable (or NaN if it does not exist)
)doc")
                        .def("get_float_list_settings",
                             &tom::ObservationAncilliarySimulationSettings::
                                 getAncilliaryDoubleVectorData,
                             py::arg("setting_type"),
                             py::arg("throw_exception") = true,
                             R"doc(
	:param setting_type:
		Type of the setting for which the value is to be returned

	:param throw_exception:
		Boolean defining whether to throw an exception if the requested setting does not exist, or does not exist as list of floating point values.

	:return:
		Value of the requested ancilliary variable (or empty list if it does not exist)
)doc");

                    m.def("doppler_ancilliary_settings",
                          &tom::getAveragedDopplerAncilliarySettings,
                          py::arg("integration_time") = 60.0, "");

                    m.def(
                        "two_way_range_ancilliary_settings",
                        &tom::getTwoWayRangeAncilliarySettings,
                        py::arg("retransmission_delay") = 0.0,
                        // py::arg("frequency_band") =
                        // tom::FrequencyBands::x_band,
                        R"doc(Factory function for creating ancilliary settings for two-way range observable.

	Factory function for creating ancilliary settings for a two-way range observable. Specifically, this
	function can be used to create settings for the retransmission delay of the observable. NOTE:
	this function is provided for convenience, and is equivalent to calling :func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_range_ancilliary_settings`
	with a single retransmission delay.


	:param retransmission_delay:
		Retransmission delay that is to be applied to the simulation of the two-way observable
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationAncilliarySimulationSettings` with the required settings.
)doc");

                    m.def(
                        "two_way_doppler_ancilliary_settings",
                        &tom::getTwoWayAveragedDopplerAncilliarySettings,
                        py::arg("integration_time") = 60.0,
                        py::arg("retransmission_delay") = 0.0,
                        R"doc(Factory function for creating ancilliary settings for two-way averaged Doppler observable.

	Factory function for creating ancilliary settings for a two-way range observable. Specifically, this
	function can be used to create settings for the retransmission delay of the observable.  NOTE:
	this function is provided for convenience, and is equivalent to calling :func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_doppler_ancilliary_settings`
	with a single retransmission delay.


	:param integration_time:
		Integration time that is to be used for the averaged Doppler observable
	:param retransmission_delay:
		Retransmission delay that is to be applied to the simulation of the two-way observable
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationAncilliarySimulationSettings` with the required settings.
)doc");

                    m.def(
                        "n_way_range_ancilliary_settings",
                        &tom::getNWayRangeAncilliarySettings,
                        py::arg("link_end_delays") = std::vector<double>(),
                        py::arg("frequency_bands") =
                            std::vector<tom::FrequencyBands>(),
                        R"doc(Factory function for creating ancilliary settings for n-way range observable.

	Factory function for creating ancilliary settings for a n-way range observable. Specifically, this
	function can be used to create settings for the retransmission delays of the observable, for each of the retransmitters.


	:param retransmission_delays:
		Retransmission delays that are to be applied to the simulation of the n-way observable. If kept empty, this results in 0 retransmission delay at each retransmitter. If defined, this list must be the same length as the number of retransmitters, and the :math:`i^{th}` entry contains the retransmission delay of the :math:`i^{th}` retrasmitter
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationAncilliarySimulationSettings` with the required settings.
)doc");

                    m.def(
                        "n_way_doppler_ancilliary_settings",
                        &tom::getNWayAveragedDopplerAncilliarySettings,
                        py::arg("integration_time") = 60.0,
                        py::arg("link_end_delays") = std::vector<double>(),
                        py::arg("frequency_bands") =
                            std::vector<tom::FrequencyBands>(),
                        R"doc(Factory function for creating ancilliary settings for n-way averaged Doppler observable.

	Factory function for creating ancilliary settings for a n-way averaged Doppler observable. Specifically, this
	function can be used to create settings for the integration time of the observable, and the  retransmission delays for each of the retransmitters.


	:param integration_time:
		Integration time that is to be used for the averaged Doppler observable
	:param retransmission_delays:
		Retransmission delays that are to be applied to the simulation of the n-way observable. If kept empty, this results in 0 retransmission delay at each retransmitter. If defined, this list must be the same length as the number of retransmitters, and the :math:`i^{th}` entry contains the retransmission delay of the :math:`i^{th}` retrasmitter
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationAncilliarySimulationSettings` with the required settings.
)doc");

                    m.def("dsn_n_way_doppler_ancilliary_settings",
                          &tom::getDsnNWayAveragedDopplerAncillarySettings,
                          py::arg("frequency_bands"),
                          py::arg("reference_frequency_band"),
                          py::arg("reference_frequency"),
                          py::arg("integration_time") = 60.0,
                          py::arg("link_end_delays") = std::vector<double>(),
                          "");

                    m.def(
                        "tabulated_simulation_settings",
                        &tss::tabulatedObservationSimulationSettings<TIME_TYPE>,
                        py::arg("observable_type"), py::arg("link_ends"),
                        py::arg("simulation_times"),
                        py::arg("reference_link_end_type") = tom::receiver,
                        py::arg("viability_settings") =
                            std::vector<std::shared_ptr<
                                tom::ObservationViabilitySettings>>(),
                        py::arg("noise_function") = nullptr,
                        py::arg("ancilliary_settings") = nullptr,
                        R"doc(Factory function for creating settings object for observation simulation, using a predefined list of observation times.

	Factory function for creating single simulation settings object, using a predefined list of observation times.
	The list of resulting observations may be reduced compared to the ``simulation_times`` provided here, as
	only observations that meet the viability settings are retained during observation simulation (these may be
	provide directly here through the ``viability_settings`` input, or added later to the resulting settings object).


	:param observable_type:
		Observable type of which observations are to be simulated.
	:param link_ends:
		Link ends for which observations are to be simulated.
	:param simulation_times:
		List of times at which to perform the observation simulation.
	:param reference_link_end_type:
		Defines the link end (via the :class:`LinkEndType`) which is used as a reference time for the observation.
	:param viability_settings:
		Settings for the creation of the viability criteria calculators, which conduct viability checks on the simulated observations.

	:param noise_function:
		Function providing the observation noise factors as a function of observation time.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.TabulatedObservationSimulationSettings` class.
)doc");


                    m.def(
                        "tabulated_simulation_settings_list",
                        &tss::createTabulatedObservationSimulationSettingsList<
                            TIME_TYPE>,
                        py::arg("link_ends_per_observable"),
                        py::arg("simulation_times"),
                        py::arg("reference_link_end_type") = tom::receiver,
                        py::arg("viability_settings") =
                            std::vector<std::shared_ptr<
                                tom::ObservationViabilitySettings>>(),
                        R"doc(Factory function for creating a list of settings object for observation simulation, using a predefined list of observation times.

	Factory function for creating multiple tabulated observation simulation settings objects in a list. This function is
	equivalent to calling the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.tabulated_simulation_settings` repeatly, with the different
	observables and link definition provided here through `link_ends_per_observable`.
	During a single call to this function, one simulation settings object is created for each combination of observable type and link geometry given by the `link_ends_per_observable` parameter.


	:param link_ends_per_observable:
		Link geometry per observable type of which observations are to be simulated.
	:param simulation_times:
		List of times at which to perform the observation simulation.
	:param reference_link_end_type:
		Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.
		The single link end specified here will be considered as the reference link end for all simulation settings object created in the function call.

	:param viability_settings:
		Settings for the creation of the viability criteria calculators, which conduct viability checks on the simulated observations.
		The single settings list given here will be considered as potential viability settings for all simulation settings object created in the function call.

	:return:
		List of :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.TabulatedObservationSimulationSettings` objects.
)doc");

                    m.def(
                        "continuous_arc_simulation_settings",
                        &tss::perArcObservationSimulationSettings<TIME_TYPE>,
                        py::arg("observable_type"), py::arg("link_ends"),
                        py::arg("start_time"), py::arg("end_time"),
                        py::arg("interval_between_observations"),
                        py::arg("arc_limiting_constraints"),
                        py::arg("minimum_arc_duration"),
                        py::arg("maximum_arc_duration"),
                        py::arg("minimum_time_between_arcs"),
                        py::arg("reference_link_end_type") = tom::receiver,
                        py::arg("additional_viability_settings") =
                            std::vector<std::shared_ptr<
                                tom::ObservationViabilitySettings>>(),
                        py::arg("noise_function") = nullptr,
                        R"doc(Factory function for creating settings object for observation simulation, using observation times according to a requirement for a continuous tracking arc.

	Factory function for creating settings object for observation simulation. Unlike the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.tabulated_simulation_settings`
	function, the resulting settings do not define the observation times explicitly. Instead, this settings object determines the observation times adaptively during the
	simulation of the observation, with the requirement that observations should be simulated over a set of contiguous arcs (if possible). The exact algorithm meets the following conditions:

	* Observations are only simulated within the time span of ``start_time`` and ``end_time``
	* A contiguous tracking arc has simulated observations separated by ``interval_between_observations``
	* Starting from ``start_time``, an observation is simulated each ``interval_between_observations``. Once an observation is unviable, as defined by
	  the ``arc_limiting_constraints`` input, it is checked whether the arc up until that point
	  is longer in duration than ``minimum_arc_duration``. If it is, the arc is added to the simulated observations. If not, the arc is discarded. In either case, a new arc is started once a
	  viable is observation is encountered
	* If the current arc reaching a duration greater than ``maximum_arc_duration``, the arc is added to the existing observations, and a new arc is started
	* If defined (e.g. if not NaN), the current observation time is incremented by ``minimum_time_between_arcs`` when an arc has been added to the observations.

	Nominally, this algorithm ensures that any arc of observations has a minimum and maximum duration. In addition, it ensures that (if desired) there is a minimum time interval
	between two tracking arcs. This behaviour can be modified by adding ``additional_viability_settings``, which are *not* used when computing the tracking arcs, but which are instead only used
	to reduce the set of simulated observations afterwards.


	:param observable_type:
		Observable type of which observations are to be simulated.
	:param link_ends:
		Link ends for which observations are to be simulated.
	:param start_time:
		First time at which an observation is to be simulated (and checked for viability).
	:param end_time:
		Maximum time at which an observation is to be simulated (and checked for viability).
	:param interval_between_observations:
		Cadence (in seconds) of subsequent observations in an arc
	:param arc_limiting_constraints:
		List of settings for the creation of the viability criteria calculators, which are used to check if an observation is viable, and define
		whether an arc should be terminated.

	:param minimum_arc_duration:
		Minimum permissible time for a tracking arc
	:param maximum_arc_duration:
		Maximum permissible time for a tracking arc
	:param minimum_time_between_arc:
		Minimum time between two tracking arcs. If NaN, this is effectively set to the ``interval_between_observations``
	:param additional_viability_settings:
		Settings for the creation of the viability criteria calculators, which conduct viability checks on the simulated observations.
		These settings are *not* used to determine whether an arc is to be terminated, but are instead applied after the arcs have been computed.

	:param noise_function:
		Function providing the observation noise factors as a function of observation time.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.TabulatedObservationSimulationSettings` class.
)doc");

                    m.def("continuous_arc_simulation_settings_list",
                          &tss::perArcObservationSimulationSettingsList<
                              TIME_TYPE>,
                          py::arg("link_ends_per_observable"),
                          py::arg("start_time"), py::arg("end_time"),
                          py::arg("interval_between_observations"),
                          py::arg("arc_limiting_constraints"),
                          py::arg("minimum_arc_duration"),
                          py::arg("maximum_arc_duration"),
                          py::arg("minimum_time_between_arcs"),
                          py::arg("reference_link_end_type") = tom::receiver,
                          py::arg("additional_viability_settings") =
                              std::vector<std::shared_ptr<
                                  tom::ObservationViabilitySettings>>(),
                          "");

                    m.def(
                        "add_noise_function_to_all",
                        py::overload_cast<
                            const std::vector<std::shared_ptr<
                                tss::ObservationSimulationSettings<
                                    TIME_TYPE>>>&,
                            const std::function<Eigen::VectorXd(const double)>>(
                            &tss::
                                addNoiseFunctionToObservationSimulationSettingsPy),
                        py::arg("observation_simulation_settings_list"),
                        py::arg("noise_amplitude"), "");

                    m.def(
                        "add_noise_function_to_observable",
                        py::overload_cast<
                            const std::vector<std::shared_ptr<
                                tss::ObservationSimulationSettings<
                                    TIME_TYPE>>>&,
                            const std::function<Eigen::VectorXd(const double)>,
                            const tom::ObservableType>(
                            &tss::
                                addNoiseFunctionToObservationSimulationSettingsPy),
                        py::arg("observation_simulation_settings_list"),
                        py::arg("noise_amplitude"), py::arg("observable_type"),
                        "");


                    m.def(
                        "add_noise_function_to_observable_for_link_ends",
                        py::overload_cast<
                            const std::vector<std::shared_ptr<
                                tss::ObservationSimulationSettings<
                                    TIME_TYPE>>>&,
                            const std::function<Eigen::VectorXd(const double)>,
                            const tom::ObservableType,
                            const tom::LinkDefinition&>(
                            &tss::
                                addNoiseFunctionToObservationSimulationSettingsPy),
                        py::arg("observation_simulation_settings_list"),
                        py::arg("noise_amplitude"), py::arg("observable_type"),
                        py::arg("link_ends"), "");


                    m.def(
                        "add_gaussian_noise_to_all",
                        py::overload_cast<
                            const std::vector<std::shared_ptr<
                                tss::ObservationSimulationSettings<
                                    TIME_TYPE>>>&,
                            const double>(
                            &tss::
                                addGaussianNoiseFunctionToObservationSimulationSettingsPy),
                        py::arg("observation_simulation_settings_list"),
                        py::arg("noise_amplitude"),
                        R"doc(Function for adding gaussian noise function to all existing observation simulation settings.

	Function for including simple time-independent and time-uncorrelated Gaussian noise function to the simulation settings of one or more observable(s).
	The noise settings are added to all :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) in the `observation_simulation_settings`
	list.

	Note: the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects are modified in-place by this function,
	and thus the function does not return anything.


	:param observation_simulation_settings:
		Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.
	:param noise_amplitude:
		Standard deviation defining the un-biased Gaussian distribution for the noise.
	:return:
		The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.

)doc");

                    m.def(
                        "add_gaussian_noise_to_observable",
                        py::overload_cast<
                            const std::vector<std::shared_ptr<
                                tss::ObservationSimulationSettings<
                                    TIME_TYPE>>>&,
                            const double, const tom::ObservableType>(
                            &tss::
                                addGaussianNoiseFunctionToObservationSimulationSettingsPy),
                        py::arg("observation_simulation_settings_list"),
                        py::arg("noise_amplitude"), py::arg("observable_type"),
                        R"doc(Function for adding gaussian noise function to existing observation simulation settings of a given observable type.

	As :func:`~tudatpy.numerical_simulation.estimation_setup.observation.add_gaussian_noise_to_all`, except that the function only adds noise to entries of the
	`observation_simulation_settings` list that matches the specified `observable_type`.


	:param observation_simulation_settings:
		Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.
	:param noise_amplitude:
		Standard deviation defining the un-biased Gaussian distribution for the noise.
	:param observable_type:
		Identifies the observable type in the observation simulation settings to which the noise is to be added.

	:return:
		The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.

)doc");


                    m.def(
                        "add_gaussian_noise_to_observable_for_link_ends",
                        py::overload_cast<
                            const std::vector<std::shared_ptr<
                                tss::ObservationSimulationSettings<
                                    TIME_TYPE>>>&,
                            const double, const tom::ObservableType,
                            const tom::LinkDefinition&>(
                            &tss::
                                addGaussianNoiseFunctionToObservationSimulationSettingsPy),
                        py::arg("observation_simulation_settings_list"),
                        py::arg("noise_amplitude"), py::arg("observable_type"),
                        py::arg("link_definition"), "");


                    m.def(
                        "add_viability_check_to_all",
                        py::overload_cast<
                            const std::vector<std::shared_ptr<
                                tss::ObservationSimulationSettings<
                                    TIME_TYPE>>>&,
                            const std::vector<std::shared_ptr<
                                tom::ObservationViabilitySettings>>&>(
                            &tss::
                                addViabilityToObservationSimulationSettingsPy),
                        py::arg("observation_simulation_settings_list"),
                        py::arg("viability_settings"),
                        R"doc(Function for including viability checks into existing observation simulation settings.

	Function for adding viability checks to the observation simulation settings, such that only observations meeting certain conditions are retained.
	The noise settings are added to all :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) in the `observation_simulation_settings`
	list.
	Note: the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects are modified in-place by this function,
	and thus the function does not return anything.


	:param observation_simulation_settings:
		Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.
	:param viability_settings:
		List of one or more :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` objects, defining the viability checks to be included.

	:return:
		The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.

)doc");


                    m.def(
                        "add_viability_check_to_observable",
                        py::overload_cast<
                            const std::vector<std::shared_ptr<
                                tss::ObservationSimulationSettings<
                                    TIME_TYPE>>>&,
                            const std::vector<std::shared_ptr<
                                tom::ObservationViabilitySettings>>&,
                            const tom::ObservableType>(
                            &tss::
                                addViabilityToObservationSimulationSettingsPy),
                        py::arg("observation_simulation_settings_list"),
                        py::arg("viability_settings"),
                        py::arg("observable_type"), "");

                    m.def(
                        "add_viability_check_to_observable_for_link_ends",
                        py::overload_cast<
                            const std::vector<std::shared_ptr<
                                tss::ObservationSimulationSettings<
                                    TIME_TYPE>>>&,
                            const std::vector<std::shared_ptr<
                                tom::ObservationViabilitySettings>>&,
                            const tom::ObservableType,
                            const tom::LinkDefinition&>(
                            &tss::
                                addViabilityToObservationSimulationSettingsPy),
                        py::arg("observation_simulation_settings_list"),
                        py::arg("viability_settings"),
                        py::arg("observable_type"), py::arg("link_ends"),
                        R"doc(Function for including viability checks into existing observation simulation settings.

	As :func:`~tudatpy.numerical_simulation.estimation_setup.observation.add_viability_check_to_all`, except that the function only adds noise to entries of the
	`observation_simulation_settings` list that matches the specified `observable_type` and `link_definition`.


	:param observation_simulation_settings:
		Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.
	:param viability_settings:
		List of one or more :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` objects, defining the viability checks to be included.

	:param observable_type:
		Identifies the observable type in the observation simulation settings for which the viability checks are to be considered.

	:param link_definition:
		Identifies the link definition in the observation simulation settings for which the viability checks are to be considered.

	:return:
		The :class

)doc");


                    m.def(
                        "add_ancilliary_settings_to_observable",
                        py::overload_cast<
                            const std::vector<std::shared_ptr<
                                tss::ObservationSimulationSettings<
                                    TIME_TYPE>>>&,
                            const std::shared_ptr<
                                tom::ObservationAncilliarySimulationSettings>&,
                            const tom::ObservableType>(
                            &tss::
                                addAncilliarySettingsToObservationSimulationSettingsPy),
                        py::arg("observation_simulation_settings_list"),
                        py::arg("ancilliary_settings"),
                        py::arg("observable_type"), "");


                    m.def(
                        "add_ancilliary_settings_to_observable_for_link_ends",
                        py::overload_cast<
                            const std::vector<std::shared_ptr<
                                tss::ObservationSimulationSettings<
                                    TIME_TYPE>>>&,
                            const std::shared_ptr<
                                tom::ObservationAncilliarySimulationSettings>&,
                            const tom::ObservableType,
                            const tom::LinkDefinition&>(
                            &tss::
                                addAncilliarySettingsToObservationSimulationSettingsPy),
                        py::arg("observation_simulation_settings_list"),
                        py::arg("ancilliary_settings"),
                        py::arg("observable_type"), py::arg("link_ends"), "");


                    py::class_<tss::ObservationDependentVariableSettings,
                               std::shared_ptr<
                                   tss::ObservationDependentVariableSettings>>(
                        m, "ObservationDependentVariableSettings",
                        R"doc(Base class for setting observation dependent variables.

	Functional (base) class for setting observation dependent variables as part of the observation output.
	Note: The associated functionality is not yet mature enough for the end user. Class is exposed for development purposes only.

)doc");

                    m.def(
                        "add_dependent_variables_to_all",
                        py::overload_cast<
                            const std::vector<std::shared_ptr<
                                tss::ObservationSimulationSettings<
                                    TIME_TYPE>>>&,
                            const std::vector<std::shared_ptr<
                                tss::ObservationDependentVariableSettings>>&,
                            const tss::SystemOfBodies&>(
                            &tss::
                                addDependentVariablesToObservationSimulationSettingsPy),
                        py::arg("observation_simulation_settings"),
                        py::arg("dependent_variable_settings"),
                        py::arg("bodies"), "");

                    m.def(
                        "add_dependent_variables_to_observable",
                        py::overload_cast<
                            const std::vector<std::shared_ptr<
                                tss::ObservationSimulationSettings<
                                    TIME_TYPE>>>&,
                            const std::vector<std::shared_ptr<
                                tss::ObservationDependentVariableSettings>>&,
                            const tss::SystemOfBodies&,
                            const tom::ObservableType>(
                            &tss::
                                addDependentVariablesToObservationSimulationSettingsPy),
                        py::arg("observation_simulation_settings"),
                        py::arg("dependent_variable_settings"),
                        py::arg("bodies"), py::arg("observable_type"), "");

                    m.def(
                        "add_dependent_variables_to_observable_for_link_ends",
                        py::overload_cast<
                            const std::vector<std::shared_ptr<
                                tss::ObservationSimulationSettings<
                                    TIME_TYPE>>>&,
                            const std::vector<std::shared_ptr<
                                tss::ObservationDependentVariableSettings>>&,
                            const tss::SystemOfBodies&,
                            const tom::ObservableType,
                            const tom::LinkDefinition&>(
                            &tss::
                                addDependentVariablesToObservationSimulationSettingsPy),
                        py::arg("observation_simulation_settings"),
                        py::arg("dependent_variable_settings"),
                        py::arg("bodies"), py::arg("observable_type"),
                        py::arg("link_ends"), "");

                    /////////////////////////////////////////////////////////////////////////////////////////////////
                    // FREQUENCIES
                    /////////////////////////////////////////////////////////////////////////////////////////////////

                    py::enum_<tom::FrequencyBands>(m, "FrequencyBands", "")
                        .value("s_band", tom::FrequencyBands::s_band)
                        .value("x_band", tom::FrequencyBands::x_band)
                        .value("ka_band", tom::FrequencyBands::ka_band)
                        .value("ku_band", tom::FrequencyBands::ku_band);

                    m.def("dsn_default_turnaround_ratios",
                          &tom::getDsnDefaultTurnaroundRatios,
                          py::arg("uplink_band"), py::arg("downlink_band"), "");

                    m.def("cassini_turnaround_ratios",
                          &tom::getCassiniTurnaroundRatio,
                          py::arg("uplink_band"), py::arg("downlink_band"), "");

                    /////////////////////////////////////////////////////////////////////////////////////////////////
                    // ODF OBSERVATIONS
                    /////////////////////////////////////////////////////////////////////////////////////////////////

                    py::class_<tom::ProcessedOdfFileContents,
                               std::shared_ptr<tom::ProcessedOdfFileContents>>(
                        m, "ProcessedOdfFileContents", "")
                        .def_property_readonly("ground_station_names",
                                               &tom::ProcessedOdfFileContents::
                                                   getGroundStationsNames,
                                               "")
                        .def_property_readonly("processed_observable_types",
                                               &tom::ProcessedOdfFileContents::
                                                   getProcessedObservableTypes,
                                               "")
                        .def_property_readonly(
                            "start_and_end_time",
                            &tom::ProcessedOdfFileContents::getStartAndEndTime,
                            "")
                        .def_property_readonly(
                            "ignored_odf_observable_types",
                            &tom::ProcessedOdfFileContents::
                                getIgnoredRawOdfObservableTypes,
                            "")
                        .def_property_readonly("ignored_ground_stations",
                                               &tom::ProcessedOdfFileContents::
                                                   getIgnoredGroundStations,
                                               "")
                        .def_property_readonly(
                            "raw_odf_data",
                            &tom::ProcessedOdfFileContents::getRawOdfData, "");

                    m.def("process_odf_data_multiple_files",
                          py::overload_cast<
                              const std::vector<std::string>&,
                              const std::string&, const bool,
                              const std::map<std::string, Eigen::Vector3d>&>(
                              &tom::processOdfData),
                          py::arg("file_names"), py::arg("spacecraft_name"),
                          py::arg("verbose") = true,
                          py::arg("earth_fixed_ground_station_positions") =
                              tss::getApproximateDsnGroundStationPositions(),
                          "");

                    m.def(
                        "process_odf_data_single_file",
                        py::overload_cast<
                            const std::string&, const std::string&, const bool,
                            const std::map<std::string, Eigen::Vector3d>&>(
                            &tom::processOdfData),
                        py::arg("file_name"), py::arg("spacecraft_name"),
                        py::arg("verbose") = true,
                        py::arg("earth_fixed_ground_station_positions") =
                            tss::getApproximateDsnGroundStationPositions(),
                        "");

                    // Create wrapper function
                    py::cpp_function getDsnDefaultTurnaroundRatios_wrapper =
                        [](tudat::observation_models::FrequencyBands band1,
                           tudat::observation_models::FrequencyBands band2) {
                            return tom::getDsnDefaultTurnaroundRatios(band1,
                                                                      band2);
                        };

                    m.def("set_odf_information_in_bodies",
                          &tom::setOdfInformationInBodies,
                          py::arg("processed_odf_file"), py::arg("bodies"),
                          py::arg("body_with_ground_stations_name") = "Earth",
                          py::arg("turnaround_ratio_function") =
                              getDsnDefaultTurnaroundRatios_wrapper,
                          "");

                    m.def(
                        "create_odf_observed_observation_collection",
                        &tom::createOdfObservedObservationCollection<double,
                                                                     TIME_TYPE>,
                        py::arg("processed_odf_file"),
                        py::arg("observable_types_to_process") =
                            std::vector<tom::ObservableType>(),
                        py::arg("start_and_end_times_to_process") =
                            std::make_pair<TIME_TYPE, TIME_TYPE>(TUDAT_NAN,
                                                                 TUDAT_NAN),
                        "");

                    //    m.def("create_odf_observation_simulation_settings_list",
                    //          &tom::createOdfObservationSimulationSettingsList<
                    //          double, TIME_TYPE >,
                    //          py::arg("observed_observation_collection"),
                    //          );

                    m.def(
                        "change_simulation_settings_observable_types",
                        &tom::
                            changeObservableTypesOfObservationSimulationSettings<
                                double, TIME_TYPE>,
                        py::arg("observation_simulation_settings"),
                        py::arg("replacement_observable_types") =
                            std::map<tom::ObservableType, tom::ObservableType>{
                                {tom::dsn_n_way_averaged_doppler,
                                 tom::n_way_differenced_range},
                                {tom::dsn_one_way_averaged_doppler,
                                 tom::one_way_differenced_range}},
                        "");

                    /////////////////////////////////////////////////////////////////////////////////////////////////
                    // Tracking Txt OBSERVATIONS
                    /////////////////////////////////////////////////////////////////////////////////////////////////

                    m.def(
                        "create_tracking_txtfile_observation_collection",
                        py::overload_cast<
                            const std::shared_ptr<
                                tudat::input_output::TrackingTxtFileContents>,
                            const std::string,
                            const std::vector<tom::ObservableType>,
                            const std::map<std::string, Eigen::Vector3d>,
                            const tom::ObservationAncilliarySimulationSettings&,
                            std::pair<TIME_TYPE, TIME_TYPE>>(
                            &tom::createTrackingTxtFileObservationCollection<
                                double, TIME_TYPE>),
                        py::arg("raw_tracking_txtfile_contents"),
                        py::arg("spacecraft_name"),
                        py::arg("observable_types_to_process") =
                            std::vector<tom::ObservableType>(),
                        py::arg("earth_fixed_ground_station_positions") =
                            tss::getApproximateDsnGroundStationPositions(),
                        py::arg("ancillary_settings") =
                            tom::ObservationAncilliarySimulationSettings(),
                        py::arg("start_and_end_times_to_process") =
                            std::make_pair<TIME_TYPE, TIME_TYPE>(TUDAT_NAN,
                                                                 TUDAT_NAN),
                        "");

                    m.def(
                        "observation_settings_from_collection",
                        py::overload_cast<std::shared_ptr<
                            tom::ObservationCollection<double, TIME_TYPE>>>(
                            &tss::
                                getObservationSimulationSettingsFromObservations<
                                    double, TIME_TYPE>),
                        py::arg("observed_observation_collection"), "");

                    //////////////////////////////////////////// DEPRECATED
                    ///////////////////////////////////////////////

                    m.def(
                        "one_way_open_loop_doppler",
                        &tom::oneWayOpenLoopDoppler, py::arg("link_ends"),
                        py::arg("light_time_correction_settings") =
                            std::vector<std::shared_ptr<
                                tom::LightTimeCorrectionSettings>>(),
                        py::arg("bias_settings") = nullptr,
                        py::arg("transmitter_proper_time_rate_settings") =
                            nullptr,
                        py::arg("receiver_proper_time_rate_settings") = nullptr,
                        py::arg("light_time_convergence_settings") =
                            std::make_shared<
                                tom::LightTimeConvergenceCriteria>(),
                        py::arg("normalized_with_speed_of_light") = false);

                    m.def("two_way_open_loop_doppler_from_one_way_links",
                          py::overload_cast<
                              const std::shared_ptr<
                                  tom::OneWayDopplerObservationSettings>,
                              const std::shared_ptr<
                                  tom::OneWayDopplerObservationSettings>,
                              const std::shared_ptr<
                                  tom::ObservationBiasSettings>>(
                              &tom::twoWayOpenLoopDoppler),
                          py::arg("uplink_doppler_settings"),
                          py::arg("downlink_doppler_settings"),
                          py::arg("bias_settings") = nullptr);

                    m.def(
                        "two_way_open_loop_doppler",
                        py::overload_cast<
                            const tom::LinkDefinition&,
                            const std::vector<std::shared_ptr<
                                tom::LightTimeCorrectionSettings>>&,
                            const std::shared_ptr<tom::ObservationBiasSettings>,
                            const std::shared_ptr<
                                tom::LightTimeConvergenceCriteria>,
                            const bool>(&tom::twoWayOpenLoopDoppler),
                        py::arg("link_ends"),
                        py::arg("light_time_correction_settings") =
                            std::vector<std::shared_ptr<
                                tom::LightTimeCorrectionSettings>>(),
                        py::arg("bias_settings") = nullptr,
                        py::arg("light_time_convergence_settings") =
                            std::make_shared<
                                tom::LightTimeConvergenceCriteria>(),
                        py::arg("normalized_with_speed_of_light") = false);

                    m.def(
                        "one_way_closed_loop_doppler",
                        py::overload_cast<
                            const tom::LinkDefinition&,
                            const std::vector<std::shared_ptr<
                                tom::LightTimeCorrectionSettings>>,
                            const std::shared_ptr<tom::ObservationBiasSettings>,
                            const std::shared_ptr<
                                tom::LightTimeConvergenceCriteria>>(
                            &tom::oneWayClosedLoopDoppler),
                        py::arg("link_ends"),
                        py::arg("light_time_correction_settings") =
                            std::vector<std::shared_ptr<
                                tom::LightTimeCorrectionSettings>>(),
                        py::arg("bias_settings") = nullptr,
                        py::arg("light_time_convergence_settings") =
                            std::make_shared<
                                tom::LightTimeConvergenceCriteria>());

                    m.def(
                        "one_way_closed_loop_doppler",
                        py::overload_cast<
                            const tom::LinkDefinition&,
                            const std::vector<std::shared_ptr<
                                tom::LightTimeCorrectionSettings>>,
                            const std::shared_ptr<tom::ObservationBiasSettings>,
                            const std::shared_ptr<
                                tom::LightTimeConvergenceCriteria>>(
                            &tom::oneWayClosedLoopDoppler),
                        py::arg("link_ends"),
                        py::arg("light_time_correction_settings") =
                            std::vector<std::shared_ptr<
                                tom::LightTimeCorrectionSettings>>(),
                        py::arg("bias_settings") = nullptr,
                        py::arg("light_time_convergence_settings") =
                            std::make_shared<
                                tom::LightTimeConvergenceCriteria>());

                    //    m.def("gaussian_noise_function",
                    //              &ts::getGaussianDistributionNoiseFunction,
                    //          py::arg("standard_deviation"),
                    //          py::arg("mean") = 0.0,
                    //          py::arg("seed") = time(NULL),
                    //          py::arg("observable_size") = 1);
                }

            }  // namespace observation
        }  // namespace estimation_setup
    }  // namespace numerical_simulation
}  // namespace tudatpy
