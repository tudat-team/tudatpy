/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#include <tudat/astro/aerodynamics/nrlmsise00Atmosphere.h>
#include <tudat/astro/aerodynamics/nrlmsise00InputFunctions.h>
#include <tudat/astro/reference_frames/referenceFrameTransformations.h>
#include <tudat/simulation/environment_setup.h>


// #include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
// #include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace trf = tudat::reference_frames;
namespace tp = tudat::physical_constants;
namespace ta = tudat::aerodynamics;
namespace tio = tudat::input_output;

namespace tudat {
    namespace simulation_setup {
        inline std::shared_ptr<AtmosphereSettings> us76AtmosphereSettings() {
            std::string atmosphereTableFile =
                paths::getAtmosphereTablesPath() +
                "/USSA1976Until100kmPer100mUntil1000kmPer1000m.dat";
            return std::make_shared<TabulatedAtmosphereSettings>(
                atmosphereTableFile);
        }
    }  // namespace simulation_setup
}  // namespace tudat
namespace tudatpy {
    namespace numerical_simulation {
        namespace environment_setup {
            namespace atmosphere {

                PYBIND11_MODULE(expose_atmosphere, m) {
                    py::module_::import(
                        "tudatpy.numerical_simulation.environment");
                    py::enum_<tss::AtmosphereDependentVariables>(
                        m, "AtmosphereDependentVariables")
                        .value("tabulated_density",
                               tss::AtmosphereDependentVariables::
                                   density_dependent_atmosphere)
                        .value("tabulated_pressure",
                               tss::AtmosphereDependentVariables::
                                   pressure_dependent_atmosphere)
                        .value("tabulated_temperature",
                               tss::AtmosphereDependentVariables::
                                   temperature_dependent_atmosphere)
                        .value("tabulated_gas_constant",
                               tss::AtmosphereDependentVariables::
                                   gas_constant_dependent_atmosphere)
                        .value("tabulated_specific_heat_ratio",
                               tss::AtmosphereDependentVariables::
                                   specific_heat_ratio_dependent_atmosphere)
                        .value("tabulated_molar_mass",
                               tss::AtmosphereDependentVariables::
                                   molar_mass_dependent_atmosphere)
                        .export_values();

                    /////////////////////////////////////////////////////////////////////////////
                    py::class_<tss::WindModelSettings,
                               std::shared_ptr<tss::WindModelSettings>>(
                        m, "WindModelSettings",
                        R"doc(Class for providing settings for wind model.

	Functional (base) class for settings of wind models that require no information in addition to their type.
	Wind model classes requiring additional information must be created using an object derived from this class.

)doc");

                    py::class_<tss::AtmosphereSettings,
                               std::shared_ptr<tss::AtmosphereSettings>>(
                        m, "AtmosphereSettings",
                        R"doc(Base class for providing settings for atmosphere model.

	Functional (base) class for settings of atmosphere models that require no information in addition to their type.
	Atmosphere model classes requiring additional information must be created using an object derived from this class.

)doc")
                        .def_property(
                            "wind_settings",
                            &tss::AtmosphereSettings::getWindSettings,
                            &tss::AtmosphereSettings::setWindSettings,
                            R"doc(Wind model settings for the atmosphere model settings object.
	)doc");

                    py::class_<
                        tss::ExponentialAtmosphereSettings,
                        std::shared_ptr<tss::ExponentialAtmosphereSettings>,
                        tss::AtmosphereSettings>(
                        m, "ExponentialAtmosphereSettings",
                        R"doc(Class for providing settings for exponential atmosphere model.

	`AtmosphereSettings` derived class for a defining the settings of an exponential atmosphere model.
)doc");

                    // unexposed this class, because there is no factory
                    // function interface yet
                    // py::class_<tss::TabulatedAtmosphereSettings,
                    //         std::shared_ptr<tss::TabulatedAtmosphereSettings>,
                    //         tss::AtmosphereSettings>(m,
                    //         "TabulatedAtmosphereSettings",


                    m.def(
                        "constant_wind_model", &tss::constantWindModelSettings,
                        py::arg("wind_velocity"),
                        py::arg("associated_reference_frame") =
                            trf::vertical_frame,
                        R"doc(Factory function for creating wind model settings with constant wind velocity.

	Factory function for settings object, defining wind model entirely from constant wind velocity in a given reference frame.


	:param wind_velocity:
		Constant wind velocity in the specified reference frame.

	:param associated_reference_frame:
		Reference frame in which constant wind velocity is defined.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.WindModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.ConstantWindModelSettings` class
)doc");

                    m.def(
                        "custom_wind_model", &tss::customWindModelSettings,
                        py::arg("wind_function"),
                        py::arg("associated_reference_frame") =
                            trf::vertical_frame,
                        R"doc(Factory function for creating wind model settings with custom wind velocity.

	Factory function for settings object, defining wind model entirely from custom wind velocity function in a given reference frame.
	The custom wind velocity has to be given as a function of altitude, longitude, latitude and time.

	.. note:: The longitude and latitude will be passed to the function in **degree** and not in radians.
	          The altitude is in meters, and the time is a Julian date in seconds since J2000.


	:param wind_velocity:
		Custom wind velocity function (w.r.t. altitude, longitude, latitude and time) in the specified reference frame.

	:param associated_reference_frame:
		Reference frame in which wind velocity is defined.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.WindModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.CustomWindModelSettings` class
)doc");


                    m.def(
                        "exponential_predefined",
                        py::overload_cast<const std::string &>(
                            &tss::exponentialAtmosphereSettings),
                        py::arg("body_name"),
                        R"doc(Factory function for creating atmospheric model settings from pre-defined exponential model.

	Factory function for settings object, defining atmosphere model from pre-defined exponential model.
	The pre-encoded properties are available for Earth and Mars, as can be seen on the table below.
	This function creates an instance of an `AtmosphereSettings` derived `ExponentialAtmosphereSettings` object.

	.. list-table:: Pre-defined exponential atmosphere model properties
	  :widths: 25 25 25 25
	  :header-rows: 1

	  * - Property
	    - Earth
	    - Mars
	    - Units
	  * - Scale Height
	    - 7.2
	    - 11.1
	    - km
	  * - Density at Zero Altitude
	    - 1.225
	    - 0.02
	    - kg/m :math:`{}^3`
	  * - Constant Temperature
	    - 246.0
	    - 215.0
	    - K
	  * - Specific Gas Constant
	    - 287.0
	    - 197.0
	    - J/kg/K
	  * - Ratio of Specific Heats
	    - 1.4
	    - 1.3
	    - --


	:param body_name:
		Body for which pre-defined model settings are to be loaded. Available bodies "Earth", "Mars".

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.ExponentialAtmosphereSettings` class
)doc");


                    m.def(
                        "exponential",
                        py::overload_cast<const double, const double,
                                          const double, const double,
                                          const double>(
                            &tss::exponentialAtmosphereSettings),
                        py::arg("scale_height"), py::arg("surface_density"),
                        py::arg("constant_temperature") = 288.15,
                        py::arg("specific_gas_constant") = tudat::
                            physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
                        py::arg("ratio_specific_heats") = 1.4,
                        R"doc(Factory function for creating atmospheric model settings from fully parametrized exponential model.

	Factory function for settings object, defining exponential atmosphere model.
	The model is solely based on an exponentially decaying density profile with a constant temperature and composition
	(i.e. independent of time, latitude and longitude).

	The user has access to a fully parametrized model, meaning that in addition to the required input parameters ``scale_height`` and ``surface_density`` (ground-level air density),
	the user can specify non-standard values for constant temperature, gas constant and specific heats ratio.


	:param scale_height:
		Scale height for density profile of atmosphere.
	:param surface_density:
		Atmospheric density at ground level.
	:param constant_temperature:
		Constant atmospheric temperature.
	:param specific_gas_constant:
		Specific gas constant for (constant) atmospheric chemical composition.
	:param ratio_specific_heats:
		Ratio of specific heats for (constant) atmospheric chemical composition.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.ExponentialAtmosphereSettings` class
)doc");

                    py::class_<ta::NRLMSISE00Atmosphere,
                               std::shared_ptr<ta::NRLMSISE00Atmosphere>>(
                        m, "NRLMSISE00Atmosphere",
                        R"doc(NRLMSISE00 atmosphere model.

                        This class uses the NRLMSISE00 model to compute the atmospheric density and temperature. The GTD7 function is used: Neutral Atmosphere Empirical Model from the surface to the lower exosphere.

                        Currently, the ideal gas law is used to compute the speed of sound and the specific heat ratio is assumed to be constant and equal to 1.4.

                        :param solar_activity_data: Solar activity data for a range of epochs as produced by tudatpy.io.read_solar_activity_data.
                        )doc")
                        .def(
                            py::init<const std::map<
                                double,
                                std::shared_ptr<
                                    tio::solar_activity::SolarActivityData>>>(),
                            py::arg("solar_activity_data"))
                        .def("get_density",
                             &ta::NRLMSISE00Atmosphere::getDensity,
                             py::arg("altitude"), py::arg("longitude"),
                             py::arg("latitude"), py::arg("time"),
                             R"doc(Get local density

                            Returns the local density at the given altitude,
                            longitude, latitude and time.

                            :param altitude: Altitude at which to get the density. [m]
                            :param longitude: Longitude at which to get the density [rad].
                            :param latitude: Latitude at which to get the density [rad].
                            :param time: Time at which density is to be computed [seconds since J2000].
                            :return: Local density. [kg/m^3]
                            )doc");

                    // POSSIBLY DEPRECATED
                    m.def(
                        "nrlmsise00", &tss::nrlmsise00AtmosphereSettings,
                        py::arg("space_weather_file") =
                            tudat::paths::getSpaceWeatherDataPath() +
                            "/sw19571001.txt",
                        R"doc(Factory function for creating NRLMSISE-00 atmospheric model settings.

	Factory function for settings object, defining atmosphere model in accordance to the NRLMSISE-00 global reference model for Earth's atmosphere.


	:param space_weather_file:
		File to be used for space weather characteristics as a function of time (e.g. F10.7, Kp, etc.). The file is typically taken from here `celestrak <https://celestrak.org/SpaceData/sw19571001.txt>`_ (note that the file in your resources path will not be the latest version of this file; download and replace your existing file if required). Documentation on the file is given `here <https://celestrak.org/SpaceData/SpaceWx-format.php>`_
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` class
)doc");

                    m.def("tabulated", &tss::tabulatedAtmosphereSettings,
                          py::arg("atmosphere_data_file"),
                          py::arg("dependent_variable_names") =
                              std::vector<ta::AtmosphereDependentVariables>(
                                  {tss::density_dependent_atmosphere,
                                   tss::pressure_dependent_atmosphere,
                                   tss::temperature_dependent_atmosphere}),
                          py::arg("specific_gas_constant") =
                              tp::SPECIFIC_GAS_CONSTANT_AIR,
                          py::arg("ratio_of_specific_heats") = 1.4);

                    m.def("us76", &tss::us76AtmosphereSettings, "");

                    m.def(
                        "custom_constant_temperature",
                        py::overload_cast<
                            const std::function<double(const double)>,
                            const double, const double, const double>(
                            &tss::customConstantTemperatureAtmosphereSettings),
                        py::arg("density_function"),
                        py::arg("constant_temperature"),
                        py::arg("specific_gas_constant") = tudat::
                            physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
                        py::arg("ratio_of_specific_heats") = 1.4,
                        R"doc(Factory function for creating atmospheric model settings from custom density profile.

	Factory function for settings object, defining constant temperature atmosphere model from custom density profile.
	The user is specifying the density profile as a function of altitude.
	The value of pressure is computed by assuming hydrostatic equilibrium, temperature, gas constant and the ratio of specific heats are modelled as constants.


	:param density_function:
		Function to retrieve the density at the current altitude.

	:param constant_temperature:
		Constant atmospheric temperature.
	:param specific_gas_constant:
		Specific gas constant for (constant) atmospheric chemical composition.
	:param ratio_specific_heats:
		Ratio of specific heats for (constant) atmospheric chemical composition.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.CustomConstantTemperatureAtmosphereSettings` class
)doc");

                    m.def(
                        "custom_four_dimensional_constant_temperature",
                        py::overload_cast<const std::function<double(
                                              const double, const double,
                                              const double, const double)>,
                                          const double, const double,
                                          const double>(
                            &tss::customConstantTemperatureAtmosphereSettings),
                        py::arg("density_function"),
                        py::arg("constant_temperature"),
                        py::arg("specific_gas_constant") = tudat::
                            physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
                        py::arg("ratio_of_specific_heats") = 1.4,
                        R"doc(Factory function for creating atmospheric model settings from custom density profile.

	Factory function for settings object, defining constant temperature atmosphere model from custom density profile.
	The user is specifying the density profile as a function of altitude, longitude, latitude and time.

	.. note:: The longitude and latitude will be passed to the function in **degree** and not in radians.
	          The altitude is in meters, and the time is a Julian date in seconds since J2000.


	:param density_function:
		Function to retrieve the density at the current altitude, longitude, latitude and time.

	:param constant_temperature:
		Constant atmospheric temperature.
	:param specific_gas_constant:
		Specific gas constant for (constant) atmospheric chemical composition.
	:param ratio_specific_heats:
		Ratio of specific heats for (constant) atmospheric chemical composition.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.CustomConstantTemperatureAtmosphereSettings` class
)doc");


                    m.def(
                        "scaled_by_function",
                        py::overload_cast<
                            const std::shared_ptr<tss::AtmosphereSettings>,
                            const std::function<double(const double)>,
                            const bool>(&tss::scaledAtmosphereSettings),
                        py::arg("unscaled_atmosphere_settings"),
                        py::arg("density_scaling_function"),
                        py::arg("is_scaling_absolute") = false,
                        R"doc(Factory function for creating scaled atmospheric model settings.

	Factory function for settings object, defining atmospheric model based on scaling an existing atmospheric settings object.
	The user can apply custom scaling factors (or absolute values) to the air densities of the existing model settings (for instance for an uncertainty analysis).


	:param unscaled_atmosphere_settings:
		Sets base settings of atmosphere model to be scaled.
	:param density_scaling_function:
		Specifies air density scaling factor as a function of time.
	:param is_scaling_absolute:
		Boolean indicating whether density scaling is absolute. Setting this boolean to true will add the scaling value to the baseline density, instead of the default behaviour of multiplying the baseline density by the scaling value.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.ScaledAtmosphereSettings` class.
)doc");

                    m.def(
                        "scaled_by_constant",
                        py::overload_cast<
                            const std::shared_ptr<tss::AtmosphereSettings>,
                            const double, const bool>(
                            &tss::scaledAtmosphereSettings),
                        py::arg("unscaled_atmosphere_settings"),
                        py::arg("density_scaling"),
                        py::arg("is_scaling_absolute") = false,
                        R"doc(Factory function for creating scaled atmospheric model settings.

	Factory function for settings object, defining atmospheric model based on an scaling of an existing atmospheric settings object.
	The user can apply a scaling factor (or an absolute value) to the air densities of the existing model settings (for instance for an uncertainty analysis).


	:param unscaled_atmosphere_settings:
		Sets base settings of atmosphere model to be scaled.
	:param density_scaling:
		Constant scaling factor to be applied to the entire air density profile.
	:param is_scaling_absolute:
		Boolean indicating whether density scaling is absolute. Setting this boolean to true will add the scaling value to the baseline density, instead of the default behaviour of multiplying the baseline density by the scaling value.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.ScaledAtmosphereSettings` class.
)doc");
                }

            }  // namespace atmosphere
        }  // namespace environment_setup
    }  // namespace numerical_simulation
}  // namespace tudatpy
