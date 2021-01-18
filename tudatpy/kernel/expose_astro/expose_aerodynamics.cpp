/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include "expose_aerodynamics.h"

#include <tudat/astro/aerodynamics/aerodynamicAcceleration.h>
#include <tudat/astro/aerodynamics/aerodynamicCoefficientGenerator.h>
#include <tudat/astro/aerodynamics/aerodynamicCoefficientInterface.h>
#include <tudat/astro/aerodynamics/aerodynamicForce.h>
#include <tudat/astro/aerodynamics/aerodynamicGuidance.h>
#include <tudat/astro/aerodynamics/aerodynamics.h>
#include <tudat/astro/aerodynamics/aerodynamicTorque.h>
#include <tudat/astro/aerodynamics/atmosphereModel.h>
#include <tudat/astro/aerodynamics/controlSurfaceAerodynamicCoefficientInterface.h>
#include <tudat/astro/aerodynamics/customAerodynamicCoefficientInterface.h>
#include <tudat/astro/aerodynamics/customConstantTemperatureAtmosphere.h>
#include <tudat/astro/aerodynamics/equilibriumWallTemperature.h>
#include <tudat/astro/aerodynamics/exponentialAtmosphere.h>
#include <tudat/astro/aerodynamics/flightConditions.h>
#include <tudat/astro/aerodynamics/hypersonicLocalInclinationAnalysis.h>
//#include "aerodynamics/nrlmsise00Atmosphere.h"
//#include "aerodynamics/nrlmsise00InputFunctions.h"
#include <tudat/astro/aerodynamics/standardAtmosphere.h>
#include <tudat/astro/aerodynamics/tabulatedAtmosphere.h>
#include <tudat/astro/aerodynamics/trimOrientation.h>
#include <tudat/astro/aerodynamics/windModel.h>

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace ta = tudat::aerodynamics;
namespace tr = tudat::reference_frames;

namespace tudat
{

namespace aerodynamics
{

class PyAerodynamicGuidance : public ta::AerodynamicGuidance {
public:
    /* Inherit the constructors */
    using AerodynamicGuidance::AerodynamicGuidance;
    using AerodynamicGuidance::currentAngleOfAttack_;
    using AerodynamicGuidance::currentAngleOfSideslip_;
    using AerodynamicGuidance::currentBankAngle_;

    void updateGuidance( const double currentTime ) override {
        PYBIND11_OVERLOAD_PURE(void, AerodynamicGuidance, updateGuidance, currentTime ); }
};

}

}

namespace tudatpy {

void expose_aerodynamics(py::module &m) {

    // TODO: Note down that Pybind11 has types registered from detatched submodules.
    //  py::class_<// Forward declaration of AerodynamicAngleCalculator from "tudat/reference_frames.h"
    //      tudat::reference_frames::AerodynamicAngleCalculator,
    //      std::shared_ptr<tudat::reference_frames::AerodynamicAngleCalculator>>
    //      AerodynamicAngleCalculator(m, "AerodynamicAngleCalculator");

    py::class_<ta::AerodynamicCoefficientInterface,
            std::shared_ptr<ta::AerodynamicCoefficientInterface>>(m, "AerodynamicCoefficientInterface", "<no_doc, only_dec>");

    py::class_<ta::AerodynamicCoefficientGenerator<3, 6>,
            std::shared_ptr<ta::AerodynamicCoefficientGenerator<3, 6>>,
            ta::AerodynamicCoefficientInterface>(m, "AerodynamicCoefficientGenerator36", "<no_doc, only_dec>");

    m.def("get_default_local_inclination_mach_points", &ta::getDefaultHypersonicLocalInclinationMachPoints,
          py::arg( "mach_regime" ) = "Full" );


    m.def("get_default_local_inclination_angle_of_attack_points", &ta::getDefaultHypersonicLocalInclinationAngleOfAttackPoints );

    m.def("get_default_local_inclination_sideslip_angle_points", &ta::getDefaultHypersonicLocalInclinationAngleOfSideslipPoints );

    py::class_<ta::HypersonicLocalInclinationAnalysis,
            std::shared_ptr<ta::HypersonicLocalInclinationAnalysis>,
            ta::AerodynamicCoefficientGenerator<3, 6>>(m, "HypersonicLocalInclinationAnalysis" )
            .def(py::init<
                 const std::vector< std::vector< double > >&,
                 const std::shared_ptr< tudat::SurfaceGeometry >,
                 const std::vector< int >&,
                 const std::vector< int >&,
                 const std::vector< bool >&,
                 const std::vector< std::vector< int > >&,
                 const double,
                 const double,
                 const Eigen::Vector3d&,
                 const bool >(),
                 py::arg("independent_variable_points"),
                 py::arg("body_shape"),
                 py::arg("number_of_lines"),
                 py::arg("number_of_points"),
                 py::arg("invert_orders"),
                 py::arg("selected_methods"),
                 py::arg("reference_area"),
                 py::arg("reference_length"),
                 py::arg("moment_reference_point"),
                 py::arg("save_pressure_coefficients") = false );

    py::class_<ta::FlightConditions,
            std::shared_ptr<ta::FlightConditions>>(m, "FlightConditions")
            .def(py::init<
                 const std::shared_ptr<tudat::basic_astrodynamics::BodyShapeModel>,
                 const std::shared_ptr<tudat::reference_frames::AerodynamicAngleCalculator>>(),
                 py::arg("shape_model"),
                 py::arg("aerodynamic_angle_calculator") = std::shared_ptr< tr::AerodynamicAngleCalculator>())
            .def("get_aerodynamic_angle_calculator", &ta::FlightConditions::getAerodynamicAngleCalculator)
            .def("update_conditions", &ta::FlightConditions::updateConditions, py::arg("current_time") )
            .def_property_readonly("aerodynamic_angle_calculator", &ta::FlightConditions::getAerodynamicAngleCalculator);

    py::class_<ta::AerodynamicGuidance, ta::PyAerodynamicGuidance,
            std::shared_ptr< ta::AerodynamicGuidance > >(m, "AerodynamicGuidance")
            .def(py::init<>())
            .def("updateGuidance", &ta::AerodynamicGuidance::updateGuidance, py::arg("current_time") )
            .def_readwrite("angle_of_attack", &ta::PyAerodynamicGuidance::currentAngleOfAttack_)
            .def_readwrite("bank_angle", &ta::PyAerodynamicGuidance::currentBankAngle_)
            .def_readwrite("sideslip_angle", &ta::PyAerodynamicGuidance::currentAngleOfSideslip_);

};

};// namespace tudatpy
