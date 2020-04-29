//
// Created by ggarrett on 23-04-20.
//

#include "trampoline_classes.hpp"
#include "expose_simulation_setup.h"
#include "docstrings.h"
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
// Ephemerides.
#include "Tudat/Mathematics/Interpolators/createInterpolator.h"
//#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"
//#include "Tudat/Astrodynamics/Ephemerides/constantEphemeris.h"

// Acceleration model types
#include "Tudat/SimulationSetup/PropagationSetup/createAccelerationModels.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

#ifndef EPHEMERIS_H
#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#define EPHEMERIS_H
#endif

#include <pybind11/pybind11.h>

// Conversion for standard types (e.g. list->vector)
#include <pybind11/stl.h>

// Limited conversion for numpy<->eigen
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

// For functional support.
#include <pybind11/functional.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace te = tudat::ephemerides;
namespace tba = tudat::basic_astrodynamics;
namespace tg = tudat::gravitation;
namespace tp = tudat::propagators;

#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/Gravitation/gravityFieldModel.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationSettings.h"

namespace toec = tudat::orbital_element_conversions;

namespace tudatpy {

    void expose_simulation_setup(py::module &m) {

        // BodySettings class
        py::class_ < tss::BodySettings,
                std::shared_ptr < tss::BodySettings > > (m, "BodySettings", tudatpy::body_settings_docstring().c_str())
                        .def_readwrite("constant_mass",
                                       &tss::BodySettings::constantMass)
                        .def_readwrite("atmosphere_settings",
                                       &tss::BodySettings::atmosphereSettings)
                        .def_readwrite("ephemeris_settings",
                                       &tss::BodySettings::ephemerisSettings)
                        .def_readwrite("gravity_field_settings",
                                       &tss::BodySettings::gravityFieldSettings)
                        .def_readwrite("rotation_model_settings",
                                       &tss::BodySettings::rotationModelSettings)
                        .def_readwrite("shape_model_settings",
                                       &tss::BodySettings::shapeModelSettings)
                        .def_readwrite("radiation_pressure_settings",
                                       &tss::BodySettings::radiationPressureSettings)
                        .def_readwrite("aerodynamic_coefficient_settings",
                                       &tss::BodySettings::aerodynamicCoefficientSettings)
                        .def_readwrite("gravity_field_variation_settings",
                                       &tss::BodySettings::gravityFieldVariationSettings)
                        .def_readwrite("ground_station_settings",
                                       &tss::BodySettings::groundStationSettings);

        // Body class
        py::class_ < tss::Body, std::shared_ptr < tss::Body > > (m, "Body")
                .def(py::init<const Eigen::Vector6d &>(),
                     py::arg("state") = Eigen::Vector6d::Zero())
                .def_property("ephemeris_frame_to_base_frame",
                              &tss::Body::getEphemerisFrameToBaseFrame,
                              &tss::Body::setEphemerisFrameToBaseFrame)
                .def_property("state",
                              &tss::Body::getState,
                              &tss::Body::setState)
                .def_property("ephemeris",
                              &tss::Body::getEphemeris,
                              &tss::Body::setEphemeris)
                .def_property("gravity_field_model",
                              &tss::Body::getGravityFieldModel,
                              &tss::Body::setGravityFieldModel);

        // getDefaultBodySettings (overload 1)
        m.def("get_default_body_settings",
              py::overload_cast<
                      const std::vector <std::string> &,
                      const double,
                      const double,
                      const double>(&tss::getDefaultBodySettings),
              py::arg("bodies"),
              py::arg("initial_time"),
              py::arg("final_time"),
              py::arg("time_step")
        );

        // getDefaultBodySettings (overload 2)
        m.def("get_default_body_settings",
              py::overload_cast<
                      const std::vector <std::string> &>(&tss::getDefaultBodySettings),
              py::arg("bodies"));

        // ephemerides (base class needs trampoline for inheritance in Python.
        py::class_ <te::Ephemeris, std::shared_ptr<te::Ephemeris>> ephemeris(m, "Ephemeris");

        // TODO: Try include base class eventually.
        py::class_ < te::ConstantEphemeris, std::shared_ptr < te::ConstantEphemeris >, te::Ephemeris > (
                m, "ConstantEphemeris")
                .def(py::init<
        const std::function<Eigen::Vector6d()>,
        const std::string &,
        const std::string &>(),
                py::arg("constant_state_function"),
                py::arg("reference_frame_origin") = "SSB",
                py::arg("reference_frame_orientation") = "ECLIPJ2000")
        .def(py::init<
                     const Eigen::Vector6d,
                     const std::string &,
                     const std::string &>(),
             py::arg("constant_state"),
             py::arg("reference_frame_origin") = "SSB",
             py::arg("reference_frame_orientation") = "ECLIPJ2000")
                .def("get_cartesian_state", &te::ConstantEphemeris::getCartesianState,
                     py::arg("seconds_since_epoch") = 0.0)
                .def("update_constant_state", &te::ConstantEphemeris::updateConstantState,
                     py::arg("new_state"))
                .def("get_cartesian_position", &te::ConstantEphemeris::getCartesianPosition,
                     py::arg("seconds_since_epoch"))
                .def("get_cartesian_velocity", &te::ConstantEphemeris::getCartesianVelocity,
                     py::arg("seconds_since_epoch"))
                .def("get_cartesian_long_state", &te::ConstantEphemeris::getCartesianLongState,
                     py::arg("seconds_since_epoch"))
                .def("get_cartesian_state_from_extended_time",
                     &te::ConstantEphemeris::getCartesianStateFromExtendedTime,
                     py::arg("current_time"))
                .def("get_cartesian_long_state_from_extended_time",
                     &te::ConstantEphemeris::getCartesianLongStateFromExtendedTime,
                     py::arg("current_time"));

        py::enum_<tss::EphemerisType>(ephemeris, "EphemerisType")
                .value("approximate_planet_positions", tss::approximate_planet_positions)
                .value("direct_spice_ephemeris", tss::direct_spice_ephemeris)
                .value("interpolated_spice", tss::interpolated_spice)
                .value("constant_ephemeris", tss::constant_ephemeris)
                .value("kepler_ephemeris", tss::kepler_ephemeris)
                .value("custom_ephemeris", tss::custom_ephemeris);


        py::class_ < tss::EphemerisSettings, std::shared_ptr < tss::EphemerisSettings >> (m, "EphemerisSettings")
                .def(py::init<const tss::EphemerisType,
                             const std::string &,
                             const std::string &>(),
                     py::arg("ephemeris_type"),
                     py::arg("frame_origin") = "SSB",
                     py::arg("frame_orientation"))
                .def("get_ephemeris_type", &tss::EphemerisSettings::getEphemerisType)
                .def("get_frame_origin", &tss::EphemerisSettings::getFrameOrigin)
                .def("get_frame_orientation", &tss::EphemerisSettings::getFrameOrientation)
                .def("get_multi_arc_ephemeris", &tss::EphemerisSettings::getMakeMultiArcEphemeris)
                .def("reset_frame_origin", &tss::EphemerisSettings::resetFrameOrigin)
                .def("reset_frame_orientation", &tss::EphemerisSettings::resetFrameOrientation)
                .def("reset_make_multi_arc_ephemeris", &tss::EphemerisSettings::resetMakeMultiArcEphemeris);

        py::class_ <
        tss::DirectSpiceEphemerisSettings,
                std::shared_ptr < tss::DirectSpiceEphemerisSettings >,
                tss::EphemerisSettings > (m, "DirectSpiceEphemerisSettings")
                        .def(py::init<const std::string,
                                     const std::string,
                                     const bool,
                                     const bool,
                                     const bool,
                                     const tss::EphemerisType>(),
                             py::arg("frame_origin") = "SSB",
                             py::arg("frame_orientation") = "ECLIPJ2000",
                             py::arg("correct_for_stellar_aberration") = false,
                             py::arg("correct_for_light_time_aberration") = false,
                             py::arg("converge_light_time_aberration") = false,
                             py::arg("ephemeris_type") = tss::direct_spice_ephemeris)
//                         py::arg("ephemeris_type") = tss::direct_spice_ephemeris)
                        .def("get_correct_for_steller_aberration",
                             &tss::DirectSpiceEphemerisSettings::getCorrectForStellarAberration)
                        .def("get_correct_for_steller_aberration",
                             &tss::DirectSpiceEphemerisSettings::getCorrectForLightTimeAberration)
                        .def("get_converge_light_time_aberration",
                                // TODO : Fix getConvergeLighTimeAberration typo in Tudat.
                             &tss::DirectSpiceEphemerisSettings::getConvergeLighTimeAberration);

        py::class_ <
        tss::InterpolatedSpiceEphemerisSettings,
                std::shared_ptr < tss::InterpolatedSpiceEphemerisSettings >,
                tss::DirectSpiceEphemerisSettings > (m, "InterpolatedSpiceEphemerisSettings")
                        .def(py::init < double,
                             double,
                             double,
                             std::string,
                             std::string,
                             std::shared_ptr < tudat::interpolators::InterpolatorSettings > > (),
                             py::arg("initial_time"),
                             py::arg("final_time"),
                             py::arg("time_step"),
                             py::arg("frame_origin") = "SSB",
                             py::arg("frame_orientation") = "ECLIPJ2000",
                             py::arg("interpolator_settings") = std::make_shared<tudat::interpolators::LagrangeInterpolatorSettings>(
                                     6));

        py::class_ <
        tss::ApproximatePlanetPositionSettings,
                std::shared_ptr < tss::ApproximatePlanetPositionSettings >,
                tss::EphemerisSettings > (m, "ApproximatePlanetPositionSettings")
                        .def(py::init<const tudat::ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData,
                                     const bool>(),
                             py::arg("body_identifier"),
                             py::arg("use_circular_coplanar_approximation"))
                        .def("get_body_identifier",
                             &tss::ApproximatePlanetPositionSettings::getBodyIdentifier)
                        .def("get_use_circular_coplanar_approximation",
                             &tss::ApproximatePlanetPositionSettings::getUseCircularCoplanarApproximation);

        py::class_ <
        tss::ConstantEphemerisSettings,
                std::shared_ptr < tss::ConstantEphemerisSettings >,
                tss::EphemerisSettings > (m, "ConstantEphemerisSettings")
                        .def(py::init<const Eigen::Vector6d &,
                                     const std::string &,
                                     const std::string &>(),
                             py::arg("constant_state"),
                             py::arg("frame_origin") = "SSB",
                             py::arg("frame_orientation") = "ECLIPJ2000");

        py::class_ <
        tss::CustomEphemerisSettings,
                std::shared_ptr < tss::CustomEphemerisSettings >,
                tss::EphemerisSettings > (m, "CustomEphemerisSettings")
                        .def(py::init<
        const std::function<Eigen::Vector6d(const double)>,
        const std::string &,
        const std::string &>(),
                py::arg("custom_state_function"),
                py::arg("frame_origin") = "SSB",
                py::arg("frame_orientation") = "ECLIPJ2000")
        .def("get_custom_state_function", &tss::CustomEphemerisSettings::getCustomStateFunction);

        py::class_ <
        tss::KeplerEphemerisSettings,
                std::shared_ptr < tss::KeplerEphemerisSettings >,
                tss::EphemerisSettings > (m, "KeplerEphemerisSettings")
                        .def(py::init<const Eigen::Vector6d &,
                                     const double,
                                     const double,
                                     const std::string &,
                                     const std::string &,
                                     const double,
                                     const double>(),
                             py::arg("initial_state_in_keplerian_elements"),
                             py::arg("epoch_of_initial_state"),
                             py::arg("central_body_gravitational_parameter"),
                             py::arg("reference_frame_origin") = "SSB",
                             py::arg("reference_frame_orientation") = "ECLIPJ2000",
                             py::arg("root_finder_absolute_tolerance") = 200.0 * std::numeric_limits<double>::epsilon(),
                             py::arg("root_finder_maximum_number_of_iterations") = 1000.0)
                        .def("get_initial_state_in_keplerian_elements",
                             &tss::KeplerEphemerisSettings::getInitialStateInKeplerianElements)
                        .def("get_epoch_of_initial_state",
                             &tss::KeplerEphemerisSettings::getEpochOfInitialState)
                        .def("get_central_body_gravitational_parameter",
                             &tss::KeplerEphemerisSettings::getCentralBodyGravitationalParameter)
                        .def("get_root_finder_absolute_tolerance",
                             &tss::KeplerEphemerisSettings::getRootFinderAbsoluteTolerance)
                        .def("get_root_finder_maximum_number_of_iterations",
                             &tss::KeplerEphemerisSettings::getRootFinderMaximumNumberOfIterations);

        py::class_ <
        tss::TabulatedEphemerisSettings,
                std::shared_ptr < tss::TabulatedEphemerisSettings >,
                tss::EphemerisSettings > (m, "TabulatedEphemerisSettings")
                        .def(py::init<const std::map<double, Eigen::Vector6d> &,
                                std::string,
                                std::string>())
                        .def("get_body_state_history",
                             &tss::TabulatedEphemerisSettings::getBodyStateHistory)
                        .def("get_use_long_double_states",
                             &tss::TabulatedEphemerisSettings::getUseLongDoubleStates)
                        .def("set_use_long_double_states",
                             &tss::TabulatedEphemerisSettings::setUseLongDoubleStates);

        // Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.cpp
        m.def("create_tabulated_ephemeris_from_spice", &tss::createTabulatedEphemerisFromSpice<>,
              py::arg("body"),
              py::arg("initial_time"),
              py::arg("end_time"),
              py::arg("time_step"),
              py::arg("observer_name"),
              py::arg("reference_frame_name"),
              py::arg("interpolator_settings") = std::make_shared<
                      tudat::interpolators::LagrangeInterpolatorSettings>(8));

        // Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.cpp
        m.def("create_body_ephemeris", &tss::createBodyEphemeris,
              py::arg("ephemeris_settings"),
              py::arg("body_name"));

        m.def("get_safe_interpolation_interval", &tss::getSafeInterpolationInterval,
              py::arg("ephemeris_model"));

        // Tudat/SimulationSetup/EnvironmentSetup/createBodies.h
        m.def("set_global_frame_body_ephemerides", &tss::setGlobalFrameBodyEphemerides < double, double > );

//            void setGlobalFrameBodyEphemerides( const NamedBodyMap& bodyMap,
//                                                const std::string& globalFrameOrigin,
//                                                const std::string& globalFrameOrientation )

        m.def("create_bodies", &tss::createBodies);

        py::class_ < tss::AccelerationSettings, std::shared_ptr <
                                                tss::AccelerationSettings >> (m, "AccelerationSettings")
                                                        .def(py::init<const tudat::basic_astrodynamics::AvailableAcceleration>(),
                                                             py::arg("acceleration_type"));



//            basic_astrodynamics::AccelerationMap createAccelerationModelsMap(
//                    const NamedBodyMap& bodyMap,
//                    const SelectedAccelerationMap& selectedAccelerationPerBody,
//                    const std::vector< std::string >& propagatedBodies,
//                    const std::vector< std::string >& centralBodies );

        // createAccelerationModelsMap (overload 1)
        m.def("create_acceleration_models_dict",
              py::overload_cast<
                      const tss::NamedBodyMap &,
                      const tss::SelectedAccelerationMap &,
                      const std::vector <std::string> &,
                      const std::vector <std::string> &
              >(&tss::createAccelerationModelsMap),
              py::arg("body_dict"),
              py::arg("selected_acceleration_per_body"),
              py::arg("propagated_bodies"),
              py::arg("central_bodies"));


//            basic_astrodynamics::AccelerationMap createAccelerationModelsMap(
//                    const NamedBodyMap& bodyMap,
//                    const SelectedAccelerationMap& selectedAccelerationPerBody,
//                    const std::map< std::string, std::string >& centralBodies );

        // createAccelerationModelsMap (overload 2)
        m.def("create_acceleration_models_dict",
              py::overload_cast<
                      const tss::NamedBodyMap &,
                      const tss::SelectedAccelerationMap &,
                      const std::map <std::string, std::string> &
              >(&tss::createAccelerationModelsMap),
              py::arg("body_dict"),
              py::arg("selected_acceleration_per_body"),
              py::arg("central_bodies"));

        // Required for python conversion.
//            TypeError: Unable to convert function return value to a Python type! The signature was
//            (body_dict: Dict[str, core.simulation_setup.Body], selected_acceleration_per_body: Dict[str, Dict[str, List[core.simulation_setup.AccelerationSettings]]], propagated_bodies: List[str], central_bodies: List[str]) -> Dict[str, Dict[str, List[tudat::basic_astrodynamics::AccelerationModel<Eigen::Matrix<double, 3, 1, 0, 3, 1> >]]]
//            py::class_<te::Ephemeris, std::shared_ptr<te::Ephemeris>> ephemeris(m, "Ephemeris");

        // TODO: Check with Dominic how many specialised AccelerationDataTypes exist.
//        py::class_ <
//        tba::AccelerationModel < Eigen::Vector3d > ,
//                std::shared_ptr < tba::AccelerationModel < Eigen::Vector3d >>
//                > acceleration_model(m, "AccelerationModel");

        // CHALLENGE
//            free(): invalid pointer
//            Process finished with exit code 134 (interrupted by signal 6: SIGABRT)
// SOLVED by including std::shared_ptr<tba::AccelerationModel<Eigen::Vector3d>

        //! Keplerian elements indices.
        //enum KeplerianElementIndices
        //{
        //    semiMajorAxisIndex = 0,
        //    eccentricityIndex = 1,
        //    inclinationIndex = 2,
        //    argumentOfPeriapsisIndex = 3,
        //    longitudeOfAscendingNodeIndex = 4,
        //    trueAnomalyIndex = 5,
        //    semiLatusRectumIndex = 0
        //};

//            namespace orbital_element_conversions {
//
//        py::enum_<toec::KeplerianElementIndices>(m, "KeplerianElementIndices")
//                .value("semi_major_axis_index", toec::KeplerianElementIndices::semiMajorAxisIndex)
//                .value("eccentricity_index", toec::KeplerianElementIndices::eccentricityIndex)
//                .value("inclination_index", toec::KeplerianElementIndices::inclinationIndex)
//                .value("argument_of_periapsis_index", toec::KeplerianElementIndices::argumentOfPeriapsisIndex)
//                .value("longitude_of_ascending_node_index",
//                       toec::KeplerianElementIndices::longitudeOfAscendingNodeIndex)
//                .value("true_anomaly_index", toec::KeplerianElementIndices::trueAnomalyIndex)
//                .value("semi_latus_rectum_index", toec::KeplerianElementIndices::semiLatusRectumIndex)
//                .export_values();

//            template< typename ScalarType = double >
//            Eigen::Matrix< ScalarType, 6, 1 > convertKeplerianToCartesianElements(
//                    const Eigen::Matrix< ScalarType, 6, 1 >& keplerianElements,
//                    const ScalarType centralBodyGravitationalParameter )

//        m.def("convert_keplerian_to_cartesian_elements", &toec::convertKeplerianToCartesianElements<>);
//            }

//        py::class_ < tg::GravityFieldModel,
//                std::shared_ptr < tg::GravityFieldModel >> (m, "GravityFieldModel")
//                        .def(py::init<
//                                     const double,
//                                     const std::function<void()>
//                             >(),
//                             py::arg("gravitational_parameter"),
//                             py::arg("update_inertia_tensor") = std::function<void()>())
//                        .def("get_gravitational_parameter", &tg::GravityFieldModel::getGravitationalParameter)
//                        .def_property("gravitational_parameter", &tg::GravityFieldModel::getGravitationalParameter,
//                                      &tg::GravityFieldModel::resetGravitationalParameter);

//
//            py::class_<tp::PropagatorSettings<double>,
//                    std::shared_ptr<tp::PropagatorSettings<double>>>(m, "PropagatorSettings")
//                    .def(py::init<
//                                 const Eigen::Matrix<double, Eigen::Dynamic, 1>,
//                                 const bool
//                         >(),
//                         py::arg("initial_body_states"),
//                         py::arg("is_multi_arc"))


//            py::class_<tp::SingleArcPropagatorSettings<double>,
//                    std::shared_ptr<tp::SingleArcPropagatorSettings<double>>>(m, "SingleArcPropagatorSettings")
//                    .def(py::init<
//                                 const tp::IntegratedStateType,
//                                 const Eigen::Matrix<double, Eigen::Dynamic, 1>,
//                                 const std::shared_ptr<tp::PropagationTerminationSettings>,
//                                 const std::shared_ptr<tp::DependentVariableSaveSettings>,
//                                 const double
//                         >(),
//                         py::arg("state_type"),
//                         py::arg("initial_body_states"),
//                         py::arg("termination_settings"),
//                         py::arg("dependent_variables_to_save") = std::shared_ptr<DependentVariableSaveSettings>(),
//                         py::arg("print_interval") = TUDAT_NAN)
//
//        py::enum_<tp::TranslationalPropagatorType>(m, "TranslationalPropagatorType")
//                .value("undefined_translational_propagator",
//                       tp::TranslationalPropagatorType::undefined_translational_propagator)
//                .value("cowell",
//                       tp::TranslationalPropagatorType::cowell)
//                .value("encke",
//                       tp::TranslationalPropagatorType::encke)
//                .value("gauss_keplerian",
//                       tp::TranslationalPropagatorType::gauss_keplerian)
//                .value("gauss_modified_equinoctial",
//                       tp::TranslationalPropagatorType::gauss_modified_equinoctial)
//                .value("unified_state_model_quaternions",
//                       tp::TranslationalPropagatorType::unified_state_model_quaternions)
//                .value("unified_state_model_modified_rodrigues_parameters",
//                       tp::TranslationalPropagatorType::unified_state_model_modified_rodrigues_parameters)
//                .value("unified_state_model_exponential_map",
//                       tp::unified_state_model_exponential_map)
//                .export_values();
//
//        py::class_ < tp::DependentVariableSaveSettings,
//                std::shared_ptr < tp::DependentVariableSaveSettings >> (m, "DependentVariableSaveSettings")
//                        .def(py::init<
//                                     const std::vector <std::shared_ptr<tp::SingleDependentVariableSaveSettings>>,
//                                     const bool
//                             >(),
//                             py::arg("dependent_variables"),
//                             py::arg("print_dependent_variable_types") = true);
//
//        py::class_ <
//        tp::PropagatorSettings < double > ,
//                std::shared_ptr < tp::PropagatorSettings < double >>
//                > PropagatorSettings_(m, "PropagatorSettings");
//
//        py::class_ <
//        tp::SingleArcPropagatorSettings < double > ,
//                std::shared_ptr < tp::SingleArcPropagatorSettings < double >>,
//                tp::PropagatorSettings < double >
//                > SingleArcPropagatorSettings_(m, "SingleArcPropagatorSettings");
//
//        py::class_ <
//        tp::TranslationalStatePropagatorSettings < double > ,
//                std::shared_ptr < tp::TranslationalStatePropagatorSettings < double >>,
//                tp::SingleArcPropagatorSettings < double >
//                > (m, "TranslationalStatePropagatorSettings")
//                        .def( // ctor 1
//                                py::init<
//                                        const std::vector <std::string> &,
//                                        const tba::AccelerationMap &,
//                                        const std::vector <std::string> &,
//                                        const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
//                                        const std::shared_ptr <tp::PropagationTerminationSettings>,
//                                        const tp::TranslationalPropagatorType,
//                                        const std::shared_ptr <tp::DependentVariableSaveSettings>,
//                                        const double>(),
//                                py::arg("central_bodies"),
//                                py::arg("accelerations_map"),
//                                py::arg("bodies_to_integrate"),
//                                py::arg("initial_body_states"),
//                                py::arg("termination_settings"),
//                                py::arg("propagator") = tp::TranslationalPropagatorType::cowell,
//                                py::arg("dependent_variables_to_save") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
//                                py::arg("print_interval") = TUDAT_NAN)
//                        .def( // ctor 2
//                                py::init<const std::vector <std::string> &,
//                                        const tss::SelectedAccelerationMap &,
//                                        const std::vector <std::string> &,
//                                        const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
//                                        const std::shared_ptr <tp::PropagationTerminationSettings>,
//                                        const tp::TranslationalPropagatorType,
//                                        const std::shared_ptr <tp::DependentVariableSaveSettings>,
//                                        const double>(),
//                                py::arg("central_bodies"),
//                                py::arg("acceleration_settings_map"),
//                                py::arg("bodies_to_integrate"),
//                                py::arg("initial_body_states"),
//                                py::arg("termination_settings"),
//                                py::arg("propagator") = tp::cowell,
//                                py::arg("dependent_variables_to_save") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
//                                py::arg("print_interval") = TUDAT_NAN)
//                        .def( // ctor 3
//                                py::init<const std::vector <std::string> &,
//                                        const tba::AccelerationMap &,
//                                        const std::vector <std::string> &,
//                                        const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
//                                        const double,
//                                        const tp::TranslationalPropagatorType,
//                                        const std::shared_ptr <tp::DependentVariableSaveSettings>,
//                                        const double>(),
//                                py::arg("central_bodies"),
//                                py::arg("accelerations_map"),
//                                py::arg("bodies_to_integrate"),
//                                py::arg("initial_body_states"),
//                                py::arg("end_time"),
//                                py::arg("propagator") = tp::cowell,
//                                py::arg("dependent_variables_to_save") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
//                                py::arg("print_interval") = TUDAT_NAN)
//                        .def( // ctor 4
//                                py::init<const std::vector <std::string> &,
//                                        const tss::SelectedAccelerationMap &,
//                                        const std::vector <std::string> &,
//                                        const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
//                                        const double,
//                                        const tp::TranslationalPropagatorType,
//                                        const std::shared_ptr <tp::DependentVariableSaveSettings>,
//                                        const double>(),
//                                py::arg("central_bodies"),
//                                py::arg("acceleration_settings_map"),
//                                py::arg("bodies_to_integrate"),
//                                py::arg("initial_body_states"),
//                                py::arg("end_time"),
//                                py::arg("propagator") = tp::cowell,
//                                py::arg("dependent_variables_to_save") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
//                                py::arg("print_interval") = TUDAT_NAN);

        // CHALLENGE
//            from core.simulation_setup import get_default_body_settings
//            ImportError: generic_type: type "TranslationalStatePropagatorSettings" referenced unknown base type "tudat::propagators::SingleArcPropagatorSettings<double>"
        // Solution: Get rid of base class reference.

        // CHALLENGE
//            ImportError: arg(): could not convert default argument into a Python object (type not registered yet?). Compile in debug mode for more information.
        // Solution 1: Add default object to exposure script.
        // Solution 2: Respect order of declaration in Python scope. (happens during import)

        // CHALLENGE
//            ImportError: /home/ggarrett/Repositories/Work/tudatBundle/build/tudatpy/src/core.so: undefined symbol: _ZN7tudatpy19basic_astrodynamics26expose_basic_astrodynamicsERN8pybind116moduleE
        // You're going to face this a lot if you are not pedantic and aware of C++
        // Solution 1: You've written a declaration which hasn't been implemented.
        // Solution 1.1: Your signature is wrong for an implemented function/class.
        // Solution 2: You haven't linked your libraries properly for your module.



    };

}