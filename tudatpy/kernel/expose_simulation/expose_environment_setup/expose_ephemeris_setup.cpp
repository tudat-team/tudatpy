/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_ephemeris_setup.h"

#include "tudatpy/docstrings.h"
#include <tudat/simulation/environment_setup.h>
#include <tudat/astro/reference_frames/referenceFrameTransformations.h>

//#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
//#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace te = tudat::ephemerides;
namespace ti = tudat::interpolators;


namespace tudatpy {
namespace simulation {
namespace environment_setup {


    void expose_ephemeris_setup(py::module &m) {

        /////////////////////////////////////////////////////////////////////////////
        // approximatePlanetPositionsBase.h
        /////////////////////////////////////////////////////////////////////////////

        py::enum_<te::ApproximatePlanetPositionsBase::BodiesWithEphemerisData>(
                m, "BodiesWithEphemerisData", "<no_doc>")
                .value("mercury", te::ApproximatePlanetPositionsBase::mercury)
                .value("venus", te::ApproximatePlanetPositionsBase::venus)
                .value("earth_moon_barycenter", te::ApproximatePlanetPositionsBase::earthMoonBarycenter)
                .value("mars", te::ApproximatePlanetPositionsBase::mars)
                .value("jupiter", te::ApproximatePlanetPositionsBase::jupiter)
                .value("saturn", te::ApproximatePlanetPositionsBase::saturn)
                .value("uranus", te::ApproximatePlanetPositionsBase::uranus)
                .value("neptune", te::ApproximatePlanetPositionsBase::neptune)
                .value("pluto", te::ApproximatePlanetPositionsBase::pluto)
                .export_values();

        /////////////////////////////////////////////////////////////////////////////
        // createEphemeris.h (complete, unverified)
        /////////////////////////////////////////////////////////////////////////////
        py::class_<tss::EphemerisSettings,
                std::shared_ptr<tss::EphemerisSettings>>(m, "EphemerisSettings",
                        get_docstring("EphemerisSettings").c_str())
//            .def(py::init<const tss::EphemerisType,
//                 const std::string &,
//                 const std::string &>(),
//                 py::arg("ephemeris_type"),
//                 py::arg("frame_origin") = "SSB",
//                 py::arg("frame_orientation") = "ECLIPJ2000")
                .def_property("frame_origin", &tss::EphemerisSettings::getFrameOrigin,
                              &tss::EphemerisSettings::resetFrameOrigin)
                .def_property("frame_orientation", &tss::EphemerisSettings::getFrameOrientation,
                              &tss::EphemerisSettings::resetFrameOrientation)
                .def_property("make_multi_arc_ephemeris", &tss::EphemerisSettings::getMakeMultiArcEphemeris,
                              &tss::EphemerisSettings::resetMakeMultiArcEphemeris)
                .def_property_readonly("ephemeris_type", &tss::EphemerisSettings::getEphemerisType);

        py::class_<tss::DirectSpiceEphemerisSettings,
                std::shared_ptr<tss::DirectSpiceEphemerisSettings>,
                tss::EphemerisSettings>(m, "DirectSpiceEphemerisSettings",
                                        get_docstring("DirectSpiceEphemerisSettings").c_str())
//           .def(py::init<const std::string, const std::string, const bool,
//                const bool, const bool, const tss::EphemerisType>(),
//                py::arg("frame_origin") = "SSB",
//                py::arg("frame_orientation") = "ECLIPJ2000",
//                py::arg("correct_for_stellar_aberration") = false,
//                py::arg("correct_for_light_time_aberration") = false,
//                py::arg("converge_light_time_aberration") = false,
//                py::arg("ephemeris_type") = tss::direct_spice_ephemeris)
                .def_property_readonly("correct_for_stellar_aberration",
                                       &tss::DirectSpiceEphemerisSettings::getCorrectForStellarAberration)
                .def_property_readonly("correct_for_light_time_aberration",
                                       &tss::DirectSpiceEphemerisSettings::getCorrectForLightTimeAberration)
                .def_property_readonly("converge_light_time_aberration",
                        // TODO : Fix getConvergeLighTimeAberration typo in Tudat.
                                       &tss::DirectSpiceEphemerisSettings::getConvergeLighTimeAberration);


        py::class_<tss::InterpolatedSpiceEphemerisSettings,
                std::shared_ptr<tss::InterpolatedSpiceEphemerisSettings>,
                tss::DirectSpiceEphemerisSettings>(m, "InterpolatedSpiceEphemerisSettings",
                                                   get_docstring("InterpolatedSpiceEphemerisSettings").c_str())
//            .def(py::init<
//                 double, double, double, std::string, std::string,
//                 std::shared_ptr<tudat::interpolators::InterpolatorSettings>>(),
//                 py::arg("initial_time"), py::arg("final_time"), py::arg("time_step"),
//                 py::arg("frame_origin") = "SSB",
//                 py::arg("frame_orientation") = "ECLIPJ2000",
//                 py::arg("interpolator_settings") = std::make_shared<
//            tudat::interpolators::LagrangeInterpolatorSettings>(6))
                .def_property_readonly("initial_time", &tss::InterpolatedSpiceEphemerisSettings::getInitialTime)
                .def_property_readonly("final_time", &tss::InterpolatedSpiceEphemerisSettings::getFinalTime)
                .def_property_readonly("time_step", &tss::InterpolatedSpiceEphemerisSettings::getTimeStep);


        py::class_<tss::ApproximatePlanetPositionSettings,
                std::shared_ptr<tss::ApproximatePlanetPositionSettings>,
                tss::EphemerisSettings>(m, "ApproximatePlanetPositionSettings",
                                        get_docstring("ApproximatePlanetPositionSettings").c_str())
//            .def(py::init<const tudat::ephemerides::ApproximatePlanetPositionsBase::
//                 BodiesWithEphemerisData,
//                 const bool>(),
//                 py::arg("body_identifier"),
//                 py::arg("use_circular_coplanar_approximation"))
                .def("get_body_identifier",
                     &tss::ApproximatePlanetPositionSettings::getBodyIdentifier)
                .def("get_use_circular_coplanar_approximation",
                     &tss::ApproximatePlanetPositionSettings::getUseCircularCoplanarApproximation);


        py::class_<tss::ConstantEphemerisSettings,
                std::shared_ptr<tss::ConstantEphemerisSettings>,
                tss::EphemerisSettings>(m, "ConstantEphemerisSettings",
                                        get_docstring("ConstantEphemerisSettings").c_str());
//            .def(py::init<const Eigen::Vector6d &,
//                 const std::string &,
//                 const std::string &>(),
//                 py::arg("constant_state"),
//                 py::arg("frame_origin") = "SSB",
//                 py::arg("frame_orientation") = "ECLIPJ2000");

        py::class_<tss::CustomEphemerisSettings,
                std::shared_ptr<tss::CustomEphemerisSettings>,
                tss::EphemerisSettings>(m, "CustomEphemerisSettings",
                                        get_docstring("CustomEphemerisSettings").c_str())
//           .def(py::init<const std::function<Eigen::Vector6d(const double)>,
//                const std::string &,
//                const std::string &>(),
//                py::arg("custom_state_function"),
//                py::arg("frame_origin") = "SSB",
//                py::arg("frame_orientation") = "ECLIPJ2000")
                .def_property_readonly("get_custom_state_function",
                     &tss::CustomEphemerisSettings::getCustomStateFunction);


        py::class_<tss::KeplerEphemerisSettings,
                std::shared_ptr<tss::KeplerEphemerisSettings>,
                tss::EphemerisSettings>(m, "KeplerEphemerisSettings",
                                        get_docstring("KeplerEphemerisSettings").c_str())
//            .def(py::init<const Eigen::Vector6d &, const double, const double,
//                 const std::string &, const std::string &, const double,
//                 const double>(),
//                 py::arg("initial_state_in_keplerian_elements"),
//                 py::arg("epoch_of_initial_state"),
//                 py::arg("central_body_gravitational_parameter"),
//                 py::arg("frame_origin") = "SSB",
//                 py::arg("frame_orientation") = "ECLIPJ2000",
//                 py::arg("root_finder_absolute_tolerance") =
//            200.0 * std::numeric_limits<double>::epsilon(),
//                 py::arg("root_finder_maximum_number_of_iterations") = 1000.0)
                .def_property_readonly("initial_state_in_keplerian_elements",
                                       &tss::KeplerEphemerisSettings::getInitialStateInKeplerianElements)
                .def_property_readonly("epoch_of_initial_state",
                                       &tss::KeplerEphemerisSettings::getEpochOfInitialState)
                .def_property_readonly("central_body_gravitational_parameter",
                                       &tss::KeplerEphemerisSettings::getCentralBodyGravitationalParameter)
                .def_property_readonly("root_finder_absolute_tolerance",
                                       &tss::KeplerEphemerisSettings::getRootFinderAbsoluteTolerance)
                .def_property_readonly("root_finder_maximum_number_of_iterations",
                                       &tss::KeplerEphemerisSettings::getRootFinderMaximumNumberOfIterations);


        py::class_<tss::TabulatedEphemerisSettings,
                std::shared_ptr<tss::TabulatedEphemerisSettings>,
                tss::EphemerisSettings>(m, "TabulatedEphemerisSettings",
                                        get_docstring("TabulatedEphemerisSettings").c_str())
//            .def(py::init<const std::map<double, Eigen::Vector6d> &, std::string,
//                 std::string>())
                .def_property_readonly("body_state_history",
                                       &tss::TabulatedEphemerisSettings::getBodyStateHistory)
                .def_property("use_long_double_states", &tss::TabulatedEphemerisSettings::getUseLongDoubleStates,
                              &tss::TabulatedEphemerisSettings::setUseLongDoubleStates);


        m.def("create_ephemeris", &tss::createBodyEphemeris,
              py::arg("ephemeris_settings"), py::arg("body_name"),
              get_docstring("create_ephemeris").c_str());


        m.def("keplerian",
              &tss::keplerEphemerisSettings,
              py::arg("initial_keplerian_state"),
              py::arg("initial_state_epoch"),
              py::arg("central_body_gravitational_parameter"),
              py::arg("frame_origin") = "SSB",
              py::arg("frame_orientation") = "ECLIPJ2000",
              py::arg("root_finder_absolute_tolerance") = 200.0 * std::numeric_limits<double>::epsilon(),
              py::arg("root_finder_maximum_iterations") = 1000.0,
              get_docstring("keplerian").c_str());

        m.def("keplerian_from_spice",
              &tss::keplerEphemerisFromSpiceSettings,
              py::arg("body"),
              py::arg("initial_state_epoch"),
              py::arg("central_body_gravitational_parameter"),
              py::arg("frame_origin") = "SSB",
              py::arg("frame_orientation") = "ECLIPJ2000",
              py::arg("root_finder_absolute_tolerance") = 200.0 * std::numeric_limits<double>::epsilon(),
              py::arg("root_finder_maximum_iterations") = 1000.0,
              get_docstring("keplerian_from_spice").c_str());


        m.def("approximate_planet_positions",
              py::overload_cast<const std::string>(&tss::approximatePlanetPositionsSettings),
              py::arg("body_name_to_use"),
              get_docstring("approximate_planet_positions").c_str());

        m.def("approximate_planet_positions",
              py::overload_cast<>(&tss::approximatePlanetPositionsSettings),
              get_docstring("approximate_planet_positions").c_str());

        m.def("direct_spice",
              py::overload_cast<const std::string, const std::string, const std::string>(
                      &tss::directSpiceEphemerisSettings),
              py::arg("frame_origin") = "SSB",
              py::arg("frame_orientation") = "ECLIPJ2000",
              py::arg("body_name_to_use") = "",
              get_docstring("direct_spice").c_str());

        m.def("interpolated_spice",
              &tss::interpolatedSpiceEphemerisSettings,
              py::arg("initial_time"),
              py::arg("final_time"),
              py::arg("time_step"),
              py::arg("frame_origin") = "SSB",
              py::arg("frame_orientation") = "ECLIPJ2000",
              py::arg("interpolator_settings") = std::make_shared<ti::LagrangeInterpolatorSettings>(6),
              py::arg("body_name_to_use") = "",
              get_docstring("interpolated_spice").c_str());

        m.def("tabulated",
              &tss::tabulatedEphemerisSettings,
              py::arg("body_state_history"),
              py::arg("frame_origin") = "SSB",
              py::arg("frame_orientation") = "ECLIPJ2000",
              get_docstring("tabulated").c_str());

        m.def("constant",
              &tss::constantEphemerisSettings,
              py::arg("constant_state"),
              py::arg("frame_origin") = "SSB",
              py::arg("frame_orientation") = "ECLIPJ2000",
              get_docstring("constant").c_str());

        m.def("scaled",
              py::overload_cast<const std::shared_ptr<tss::EphemerisSettings>,
                      const double, const bool>(&tss::scaledEphemerisSettings),
              py::arg("unscaled_ephemeris_settings"),
              py::arg("scaling_constant"),
              py::arg("is_scaling_absolute") = false,
              get_docstring("scaled").c_str());

        m.def("scaled",
              py::overload_cast<const std::shared_ptr<tss::EphemerisSettings>,
                      const Eigen::Vector6d, const bool>(&tss::scaledEphemerisSettings),
              py::arg("unscaled_ephemeris_settings"),
              py::arg("scaling_vector"),
              py::arg("is_scaling_absolute") = false,
              get_docstring("scaled").c_str());

        m.def("scaled",
              py::overload_cast<const std::shared_ptr<tss::EphemerisSettings>,
                      const std::function<Eigen::Vector6d(const double)>, const bool>(&tss::scaledEphemerisSettings),
              py::arg("unscaled_ephemeris_settings"),
              py::arg("scaling_vector_function"),
              py::arg("is_scaling_absolute") = false,
              get_docstring("scaled").c_str());

        m.def("custom",
              &tss::customEphemerisSettings,
              py::arg("custom_state_function"),
              py::arg("frame_origin") = "SSB",
              py::arg("frame_orientation") = "ECLIPJ2000",
              get_docstring("custom").c_str());
    }



}// namespace environment_setup
}// namespace simulation
}// namespace tudatpy