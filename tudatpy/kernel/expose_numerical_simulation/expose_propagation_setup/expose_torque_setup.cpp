/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_torque_setup.h"
#include <tudat/basics/deprecationWarnings.h>

#include "tudatpy/docstrings.h"
#include <tudat/simulation/propagation_setup.h>

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

namespace tudat
{
namespace simulation_setup
{
inline std::shared_ptr< TorqueSettings > customTorqueSettingsDeprecated(
            const std::function< Eigen::Vector3d( const double ) > torqueFunction,
        const std::function< double( const double ) > scalingFunction = nullptr )
{
    static bool isWarningPrinted = false;
    if( isWarningPrinted == false )
    {
        tudat::utilities::printDeprecationWarning( "tudatpy.numerical_simulation.propagation_setup.acceleration.custom",
                             "tudatpy.numerical_simulation.propagation_setup.acceleration.custom_torque");
        isWarningPrinted = true;
    }

    return customTorqueSettings( torqueFunction, scalingFunction     );
}
}
}
namespace tudatpy {
namespace numerical_simulation {
namespace propagation_setup {
namespace torque {

    void expose_torque_setup(py::module &m) {

        py::enum_<tba::AvailableTorque>(m, "AvailableTorque",
                                        get_docstring("AvailableTorque").c_str())
                .value("torque_free_type", tba::AvailableTorque::torque_free,
                       get_docstring("AvailableTorque.torque_free_type").c_str())
                .value("underfined_type", tba::AvailableTorque::underfined_torque,
                       get_docstring("AvailableTorque.underfined_type").c_str())
                .value("second_order_gravitational_type", tba::AvailableTorque::second_order_gravitational_torque,
                       get_docstring("AvailableTorque.second_order_gravitational_type").c_str())
                .value("aerodynamic_type", tba::AvailableTorque::aerodynamic_torque,
                       get_docstring("AvailableTorque.aerodynamic_type").c_str())
                .value("spherical_harmonic_gravitational_type",
                       tba::AvailableTorque::spherical_harmonic_gravitational_torque,
                       get_docstring("AvailableTorque.spherical_harmonic_gravitational_type").c_str())
                .value("inertial_type", tba::AvailableTorque::inertial_torque,
                       get_docstring("AvailableTorque.inertial_type").c_str())
                .value("dissipative_type", tba::AvailableTorque::dissipative_torque,
                       get_docstring("AvailableTorque.dissipative_type").c_str())
                .export_values();

        py::class_<tss::TorqueSettings,
                std::shared_ptr<tss::TorqueSettings>>(m, "TorqueSettings",
                                                      get_docstring("TorqueSettings").c_str());

        py::class_<tss::SphericalHarmonicTorqueSettings,
                std::shared_ptr<tss::SphericalHarmonicTorqueSettings>,
                tss::TorqueSettings>(m, "SphericalHarmonicTorqueSettings",
                                     get_docstring("SphericalHarmonicTorqueSettings").c_str());


        m.def("aerodynamic", &tss::aerodynamicTorque,
              get_docstring("aerodynamic").c_str());


        m.def("second_degree_gravitational", &tss::secondDegreeGravitationalTorque,
              get_docstring("second_degree_gravitational").c_str());


        m.def("spherical_harmonic_gravitational", &tss::sphericalHarmonicGravitationalTorque,
              py::arg("maximum_degree"),
              py::arg("maximum_order"),
              get_docstring("spherical_harmonic_gravitational").c_str());



        m.def("custom_torque", &tss::customTorqueSettings,
                  py::arg("torque_function"),
                  py::arg("scaling_function") = nullptr,
              get_docstring("custom").c_str());

        m.def("custom", &tss::customTorqueSettingsDeprecated,
              py::arg("torque_function"),
              py::arg("scaling_function") = nullptr );

        // NOTE: the only unexposed torque model is dissipativeTorque, but it is probably obsolete
    }

}// namespace torque
}// namespace propagation_setup
}// namespace numerical_simulation
}// namespace tudatpy
