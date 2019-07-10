// TODO: Disable warnings according to detected compiler. GCC/clang.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#pragma GCC diagnostic ignored "-Wpragmas"

#include <vector>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/numpy.hpp>
#include <boost/python/list.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/stl_iterator.hpp>

#include <iostream>
#include <vector>
#include <Eigen/Core>
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"

// Ephemerides.
#include "Tudat/Mathematics/Interpolators/createInterpolator.h"
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"
#include "Tudat/Astrodynamics/Ephemerides/constantEphemeris.h"
#ifndef EPHEMERIS_H
#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#define EPHEMERIS_H

#include "conversion.h"
#include "eigen_numpy.h"
#include "typedefs.h"

#include "boost/shared_ptr.hpp"
#include "boost/python.hpp"
#include "boost/python/stl_iterator.hpp"

namespace bp = boost::python;
namespace bpn = boost::python::numpy;
namespace tss = tudat::simulation_setup;
namespace te = tudat::ephemerides;
namespace ti = tudat::interpolators;

//typedef double TimeType;
typedef Eigen::Matrix<double, 6, 1> StateType;

//! Typedef for state interpolator
typedef std::shared_ptr <ti::OneDimensionalInterpolator
        <double, StateType>> StateInterpolatorPointer;

std::map <std::string, std::shared_ptr<tss::BodySettings>>
(*getDefaultBodySettings1)(const std::vector <std::string> &, const double, const double, const double)
= &tss::getDefaultBodySettings;

std::map <std::string, std::shared_ptr<tss::BodySettings>>
(*getDefaultBodySettings2)(const std::vector <std::string> &)
= &tss::getDefaultBodySettings;

std::map <std::string, std::shared_ptr<tss::BodySettings>>
getDefaultBodySettings1Wrapped(bp::list &body_list, double initial_time, double final_time, double time_step) {
    return getDefaultBodySettings1(tudatpy::conversion::py_list_to_std_vector<std::string>(body_list),
                                   initial_time, final_time, time_step);
}

std::map <std::string, std::shared_ptr<tss::BodySettings>> getDefaultBodySettings2Wrapped(bp::list &body_list) {
    return getDefaultBodySettings2(tudatpy::conversion::py_list_to_std_vector<std::string>(body_list));
}

tss::NamedBodyMap createBodiesWrapped(stdBodySettingsPointerMap bodySettings) {
    return tss::createBodies(bodySettings);
}


void register_environment_converters() {
    // register the bodysettings-map-to-python-dict.
    boost::python::to_python_converter < std::map < std::string, stdBodySettingsPointer >,
            map_to_python_dict < std::string, stdBodySettingsPointer >> ();

    // register the dict-to-bodysettings-map.
    from_python_dict < std::map < std::string, stdBodySettingsPointer >> ();

    // register the unordered-body-map-to-python-dict converter.
    boost::python::to_python_converter <
    std::unordered_map < std::string, std::shared_ptr < tss::Body >>,
            unordered_map_to_python_dict < std::string, std::shared_ptr < tss::Body >> > ();

    // register to python-dict-to-undordered-map converter.
    from_python_dict < std::unordered_map < std::string, std::shared_ptr < tss::Body >> > ();
};


//class EphemerisWrap : Ephemeris, wrapper<Ephemeris> {
//
//
//};


BOOST_PYTHON_MODULE (_environment) {
        SetupEigenConverters();
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        //  ███████╗███╗   ██╗██╗   ██╗██╗██████╗  ██████╗ ███╗   ██╗███╗   ███╗███████╗███╗   ██╗████████╗
        //  ██╔════╝████╗  ██║██║   ██║██║██╔══██╗██╔═══██╗████╗  ██║████╗ ████║██╔════╝████╗  ██║╚══██╔══╝
        //  █████╗  ██╔██╗ ██║██║   ██║██║██████╔╝██║   ██║██╔██╗ ██║██╔████╔██║█████╗  ██╔██╗ ██║   ██║
        //  ██╔══╝  ██║╚██╗██║╚██╗ ██╔╝██║██╔══██╗██║   ██║██║╚██╗██║██║╚██╔╝██║██╔══╝  ██║╚██╗██║   ██║
        //  ███████╗██║ ╚████║ ╚████╔╝ ██║██║  ██║╚██████╔╝██║ ╚████║██║ ╚═╝ ██║███████╗██║ ╚████║   ██║
        //  ╚══════╝╚═╝  ╚═══╝  ╚═══╝  ╚═╝╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚═╝     ╚═╝╚══════╝╚═╝  ╚═══╝   ╚═╝
        //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        Py_Initialize();
        /**
         * TODO: Find out why the following desired design:
         *      from tudatpy.spice import load_standard_spice_kernels()
         *      does not correctly call tudat::spice_interface::loadStandardSpiceKernels() in the Python environment.
         */
        tudat::spice_interface::loadStandardSpiceKernels();

        /**
         * Register type converters between C++ and Python environment. ~ tudatpy/include/conversion.h
         */
        register_environment_converters();

        // Documentation for environment module.
        bp::scope().attr("__doc__") = "ENVIRONMENT"
        "Models for the physical environment are one of the cornerstones of a numerical"
        "astrodynamics toolbox. Here, we define the environment in the broadest sense, "
        "including all physical properties of the solar system, such as atmospheres and "
        "gravity fields, but also any models for the orbits and rotations of these bodies."
        "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html";

        // Expose tss::BodySettings, managed as std::shared_ptr<tss::BodySettings> in Python environment.
        bp::class_<tss::BodySettings, std::shared_ptr<tss::BodySettings> >("BodySettings")
        .add_property("constant_mass", &tss::BodySettings::constantMass, "Constant mass.")
        .add_property("atmosphere_settings", &tss::BodySettings::atmosphereSettings,
        "Settings for the atmosphere model that the body is to contain.")
        .add_property("ephemeris_settings", &tss::BodySettings::ephemerisSettings,
        "Settings for the ephemeris model that the body is to contain.")
        .add_property("gravity_field_settings", &tss::BodySettings::gravityFieldSettings,
        "Settings for the gravity field model that the body is to contain.")
        .add_property("rotation_model_settings", &tss::BodySettings::rotationModelSettings,
        "Settings for the rotation model that the body is to contain.")
        .add_property("shape_model_settings", &tss::BodySettings::shapeModelSettings,
        "Settings for the shape model that the body is to contain.")
        .add_property("radiation_pressure_settings", &tss::BodySettings::radiationPressureSettings,
        "Settings for the radiations pressure interfaces that the body is to contain (source body as key).")
        .add_property("aerodynamic_coefficient_settings", &tss::BodySettings::aerodynamicCoefficientSettings,
        "Settings for the aerodynamic coefficients that the body is to contain.")
        .add_property("gravity_field_variation_settings", &tss::BodySettings::gravityFieldVariationSettings,
        "Settings for variations of the gravity field of the body.")
        .add_property("ground_station_settings", &tss::BodySettings::groundStationSettings,
        "Settings for ground station information of the body.");

        // Expose EphemerisType enumeration to Python environment.
        bp::enum_<tss::EphemerisType>("EphemerisType")
        .value("approximate_planet_positions", tss::approximate_planet_positions)
        .value("direct_spice_ephemeris", tss::direct_spice_ephemeris)
        .value("interpolated_spice", tss::interpolated_spice)
        .value("constant_ephemeris", tss::constant_ephemeris)
        .value("kepler_ephemeris", tss::kepler_ephemeris)
        .value("custom_ephemeris", tss::custom_ephemeris);

        // Expose tss::EphemerisSettings, managed as std::shared_ptr<tss::EphemerisSettings> in Python environment.
        bp::class_<tss::EphemerisSettings, std::shared_ptr<tss::EphemerisSettings> >("EphemerisSettings",
        bp::init<tss::EphemerisType const, std::string const &, std::string const &>())
        .add_property("ephemeris_type",
                &tss::EphemerisSettings::getEphemerisType,
                "(EphemerisType) Ephemeris type.")
        .add_property("frame_origin",
                &tss::EphemerisSettings::getFrameOrigin,
                "(str) Frame origin.")
        .add_property("frame_orientation",
                &tss::EphemerisSettings::getFrameOrientation,
                "(str) Frame orientation.")
        .add_property("multi_arc_ephemeris",
                &tss::EphemerisSettings::getMakeMultiArcEphemeris,
                "(Bool) Multi-arc ephemeris.")
        .def("reset_frame_origin",
                &tss::EphemerisSettings::resetFrameOrigin,
                "f(string frame_origin) Reset frame origin.")
        .def("reset_frame_orientation",
                &tss::EphemerisSettings::resetFrameOrientation,
                "f(string frame_orientation) Reset frame orientation.")
        .def("reset_make_multi_arc_ephemeris",
                &tss::EphemerisSettings::resetMakeMultiArcEphemeris,
                "f(bool multi_arc_ephemeris) Reset make multi arc ephemeris.");

        // ["Tudat/Astrodynamics/Ephemerides/constantEphemeris.h"] Expose tss::ConstantEphemerisSettings.
        bp::class_<tss::ConstantEphemerisSettings,
                bp::bases<tss::EphemerisSettings>,
                std::shared_ptr<tss::ConstantEphemerisSettings >>(
        "ConstantEphemerisSettings", bp::init<const Eigen::Vector6d&, const std::string& , const std::string&>())
        .add_property("get_constant_state",
                &tss::ConstantEphemerisSettings::getConstantState)

        // Expose tss::CustomEphemerisSettings, managed as std::shared_ptr<tss::CustomEphemerisSettings> in Python environment.
//        bp::class_<
//        tss::CustomEphemerisSettings,
//        bp::bases<tss::EphemerisSettings>,
//        std::shared_ptr<tss::CustomEphemerisSettings >>(
//        "CustomEphemerisSettings", bp::init<const std::function< Eigen::Vector6d( const double ) > , const std::string&, const std::string&>())
//        .add_property("get_custom_state_function", &tss::CustomEphemerisSettings::getConstantState)


        // ["Tudat/SimulationSetup/EnvironmentSetup/body.h"] Expose tss::Body.
        bp::class_<tss::Body, std::shared_ptr<tss::Body> >("Body", bp::init<>())
        .def(bp::init<const Eigen::Vector6d&>())
        .add_property("ephemeris_frame_to_base_frame",
                &tss::Body::getEphemerisFrameToBaseFrame,
                &tss::Body::setEphemerisFrameToBaseFrame)
        .add_property("state",
                &tss::Body::getState,
                &tss::Body::setState)
        .add_property("ephemeris",
                &tss::Body::getEphemeris,
                &tss::Body::setEphemeris)
        ;


        bp::def("get_default_body_settings", getDefaultBodySettings1Wrapped);

        // [Tudat/SimulationSetup/EnvironmentSetup/createBodies.h]
        bp::def("set_global_frame_body_ephemerides", tss::setGlobalFrameBodyEphemerides<double, double>);

        bp::def("get_default_body_settings", getDefaultBodySettings2Wrapped);

        bp::def("create_bodies", createBodiesWrapped);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        //  ██████╗ ██████╗  ██████╗ ██████╗  █████╗  ██████╗  █████╗ ████████╗██╗ ██████╗ ███╗   ██╗
        //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔══██╗██╔════╝ ██╔══██╗╚══██╔══╝██║██╔═══██╗████╗  ██║
        //  ██████╔╝██████╔╝██║   ██║██████╔╝███████║██║  ███╗███████║   ██║   ██║██║   ██║██╔██╗ ██║
        //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══██║██║   ██║██╔══██║   ██║   ██║██║   ██║██║╚██╗██║
        //  ██║     ██║  ██║╚██████╔╝██║     ██║  ██║╚██████╔╝██║  ██║   ██║   ██║╚██████╔╝██║ ╚████║
        //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚═╝ ╚═════╝ ╚═╝  ╚═══╝
        //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // [Tudat/SimulationSetup/PropagationSetup/accelerationSettings.h]
        bp::class_<tss::AccelerationSettings,
        std::shared_ptr<tss::AccelerationSettings >>(
        "AccelerationSettings", bp::init<const Eigen::Vector6d&, const std::string& , const std::string&>())
        .add_property("get_constant_state",
        &tss::ConstantEphemerisSettings::getConstantState)



        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        //   █████╗ ███████╗████████╗██████╗  ██████╗ ██████╗ ██╗   ██╗███╗   ██╗ █████╗ ███╗   ███╗██╗ ██████╗███████╗
        //  ██╔══██╗██╔════╝╚══██╔══╝██╔══██╗██╔═══██╗██╔══██╗╚██╗ ██╔╝████╗  ██║██╔══██╗████╗ ████║██║██╔════╝██╔════╝
        //  ███████║███████╗   ██║   ██████╔╝██║   ██║██║  ██║ ╚████╔╝ ██╔██╗ ██║███████║██╔████╔██║██║██║     ███████╗
        //  ██╔══██║╚════██║   ██║   ██╔══██╗██║   ██║██║  ██║  ╚██╔╝  ██║╚██╗██║██╔══██║██║╚██╔╝██║██║██║     ╚════██║
        //  ██║  ██║███████║   ██║   ██║  ██║╚██████╔╝██████╔╝   ██║   ██║ ╚████║██║  ██║██║ ╚═╝ ██║██║╚██████╗███████║
        //  ╚═╝  ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝    ╚═╝   ╚═╝  ╚═══╝╚═╝  ╚═╝╚═╝     ╚═╝╚═╝ ╚═════╝╚══════╝
        //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // tudatBundle/tudat/Tudat/Astrodynamics/Ephemerides
        // Base class to be inherited.
        bp::class_< te::Ephemeris, boost::noncopyable, std::shared_ptr<te::Ephemeris >>("Ephemeris", bp::no_init);


        bp::class_< te::TabulatedCartesianEphemeris<double, double>, bp::bases<te::Ephemeris>, std::shared_ptr<te::TabulatedCartesianEphemeris<double, double> > >(
            "TabulatedCartesianEphemeris", bp::init<const StateInterpolatorPointer, const std::string, const std::string>())
            .def("get_cartesian_state", &te::TabulatedCartesianEphemeris<double, double>::getCartesianState,(bp::arg("secondsSinceEpoch")))
    //        .property("ephemeris", &tss::Body::getEphemeris, &tss::Body::setEphemeris)
            ;




//        "(EphemerisType) Ephemeris type.")
//        .add_property("frame_origin", &tss::EphemerisSettings::getFrameOrigin,
//        "(str) Frame origin.")
//        .add_property("frame_orientation", &tss::EphemerisSettings::getFrameOrientation,
//        "(str) Frame orientation.")
//        .add_property("multi_arc_ephemeris", &tss::EphemerisSettings::getMakeMultiArcEphemeris,
//        "(Bool) Multi-arc ephemeris.")
//        .def("reset_frame_origin", &tss::EphemerisSettings::resetFrameOrigin,
//        "f(string frame_origin) Reset frame origin.")
//        .def("reset_frame_orientation", &tss::EphemerisSettings::resetFrameOrientation,
//        "f(string frame_orientation) Reset frame orientation.")
//        .def("reset_make_multi_arc_ephemeris", &tss::EphemerisSettings::resetMakeMultiArcEphemeris,
//        "f(bool multi_arc_ephemeris) Reset make multi arc ephemeris.")



//    typedef std::unordered_map<std::string, std::shared_ptr<tss::Body> > NamedBodyMap;

//    bp::class_<tss::BodySettings, boostBodySettingsPointer>("BodySettings");

//    bp::class_<boostBodySettingsPointerMap>("DefaultBodySettings", bp::no_init)
//            .def("__init__", bp::make_constructor(&PythonWrapper::DefaultBodyMapSettings1))
//            .def("__init__", bp::make_constructor(&PythonWrapper::DefaultBodyMapSettings2));



        //                bp::class_<tss::BodySettings, bp::shared_ptr<tss::BodySettings>>("BodySettings", bp::init<Eigen::Vector6d const &>())
//                .def("__init__", Ei)
//                .def("get_ephemeris_frame_to_base_frame", &tss::Body::getEphemerisFrameToBaseFrame)
//                ;



}

#pragma GCC diagnostic pop