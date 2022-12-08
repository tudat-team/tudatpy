/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_atmosphere_setup.h"

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
namespace trf = tudat::reference_frames;
namespace tp = tudat::physical_constants;


namespace tudatpy {
namespace numerical_simulation {
namespace environment_setup {
namespace atmosphere {

    void expose_atmosphere_setup(py::module &m) {



        py::enum_<tss::AtmosphereDependentVariables>(m, "AtmosphereDependentVariables" )
                .value("tabulated_density", tss::AtmosphereDependentVariables::density_dependent_atmosphere )
                .value("tabulated_pressure", tss::AtmosphereDependentVariables::pressure_dependent_atmosphere )
                .value("tabulated_temperature", tss::AtmosphereDependentVariables::temperature_dependent_atmosphere )
                .value("tabulated_gas_constant", tss::AtmosphereDependentVariables::gas_constant_dependent_atmosphere )
                .value("tabulated_specific_heat_ratio", tss::AtmosphereDependentVariables::specific_heat_ratio_dependent_atmosphere )
                .value("tabulated_molar_mass", tss::AtmosphereDependentVariables::molar_mass_dependent_atmosphere )
                .export_values();

        /////////////////////////////////////////////////////////////////////////////
        py::class_<tss::WindModelSettings,
                std::shared_ptr<tss::WindModelSettings>>(m, "WindModelSettings",
                                                         get_docstring("WindModelSettings").c_str());

        py::class_<tss::AtmosphereSettings,
                std::shared_ptr<tss::AtmosphereSettings>>(m, "AtmosphereSettings",
                                                          get_docstring("AtmosphereSettings").c_str())
                .def_property("wind_settings", &tss::AtmosphereSettings::getWindSettings,
                              &tss::AtmosphereSettings::setWindSettings,
                              get_docstring("AtmosphereSettings.wind_settings").c_str());

        py::class_<tss::ExponentialAtmosphereSettings,
                std::shared_ptr<tss::ExponentialAtmosphereSettings>,
                tss::AtmosphereSettings>(m, "ExponentialAtmosphereSettings",
                                         get_docstring("ExponentialAtmosphereSettings").c_str());

        // unexposed this class, because there is no factory function interface yet
        // py::class_<tss::TabulatedAtmosphereSettings,
        //         std::shared_ptr<tss::TabulatedAtmosphereSettings>,
        //         tss::AtmosphereSettings>(m, "TabulatedAtmosphereSettings",
        //                                  get_docstring("TabulatedAtmosphereSettings").c_str());


        m.def("constant_wind_model",
              &tss::constantWindModelSettings,
              py::arg("wind_velocity"),
              py::arg("associated_reference_frame") = trf::vertical_frame,
              get_docstring("constant_wind_model").c_str());

        m.def("custom_wind_model",
              &tss::customWindModelSettings,
              py::arg("wind_function"),
              py::arg("associated_reference_frame") = trf::vertical_frame,
              get_docstring("custom_wind_model").c_str());


        m.def("exponential_predefined",
              py::overload_cast<const std::string &>(
                      &tss::exponentialAtmosphereSettings),
              py::arg("body_name"),
              get_docstring("exponential_predefined").c_str());


        m.def("exponential",
              py::overload_cast<const double, const double, const double,
                      const double, const double>(&tss::exponentialAtmosphereSettings),
              py::arg("scale_height"),
              py::arg("surface_density"),
              py::arg("constant_temperature") = 288.15,
              py::arg("specific_gas_constant") = tudat::physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
              py::arg("ratio_specific_heats") = 1.4,
              get_docstring("exponential").c_str());


        m.def("nrlmsise00",
              &tss::nrlmsise00AtmosphereSettings,
              py::arg("space_weather_file" ) = tudat::paths::getSpaceWeatherDataPath( ) + "/sw19571001.txt",
              get_docstring("nrlmsise00").c_str());

        m.def("tabulated",
              &tss::tabulatedAtmosphereSettings,
              py::arg("atmosphere_data_file" ),
              py::arg("dependent_variable_names" ) =
        { tss::density_dependent_atmosphere, tss::pressure_dependent_atmosphere, tss::temperature_dependent_atmosphere },
              py::arg("specific_gas_constant" ) = tp::SPECIFIC_GAS_CONSTANT_AIR,
              py::arg("ratio_of_specific_heats" ) = 1.4 );

        m.def("custom_constant_temperature",
              py::overload_cast<const std::function<double(const double)>,
                      const double, const double, const double>(&tss::customConstantTemperatureAtmosphereSettings),
              py::arg("density_function"),
              py::arg("constant_temperature"),
              py::arg("specific_gas_constant") = tudat::physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
              py::arg("ratio_of_specific_heats") = 1.4,
              get_docstring("custom_constant_temperature").c_str());

        m.def("custom_four_dimensional_constant_temperature",
              py::overload_cast<const std::function<double(const double, const double, const double, const double)>,
                      const double, const double, const double>(&tss::customConstantTemperatureAtmosphereSettings),
              py::arg("density_function"),
              py::arg("constant_temperature"),
              py::arg("specific_gas_constant") = tudat::physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
              py::arg("ratio_of_specific_heats") = 1.4,
              get_docstring("custom_four_dimensional_constant_temperature").c_str());


        m.def("scaled_by_function",
              py::overload_cast<const std::shared_ptr<tss::AtmosphereSettings>,
                      const std::function<double(const double)>, const bool>(&tss::scaledAtmosphereSettings),
              py::arg("unscaled_atmosphere_settings"),
              py::arg("density_scaling_function"),
              py::arg("is_scaling_absolute") = false,
              get_docstring("scaled_by_function").c_str());

        m.def("scaled_by_constant",
              py::overload_cast<const std::shared_ptr<tss::AtmosphereSettings>,
                      const double, const bool>(&tss::scaledAtmosphereSettings),
              py::arg("unscaled_atmosphere_settings"),
              py::arg("density_scaling"),
              py::arg("is_scaling_absolute") = false,
              get_docstring("scaled_by_constant").c_str());


    }

}// namespace atmosphere
}// namespace environment_setup
}// namespace numerical_simulation
}// namespace tudatpy
