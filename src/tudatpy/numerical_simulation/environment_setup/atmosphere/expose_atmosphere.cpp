/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_atmosphere.h"

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

namespace tudat
{
namespace simulation_setup
{
inline std::shared_ptr< AtmosphereSettings > us76AtmosphereSettings( )
{
    std::string atmosphereTableFile =
            paths::getAtmosphereTablesPath( ) + "/USSA1976Until100kmPer100mUntil1000kmPer1000m.dat";
    return std::make_shared< TabulatedAtmosphereSettings >( atmosphereTableFile );
}
}  // namespace simulation_setup
}  // namespace tudat
namespace tudatpy
{
namespace numerical_simulation
{
namespace environment_setup
{
namespace atmosphere
{

void expose_atmosphere_setup( py::module &m )
{
    // NRLMSISE00
    py::class_< ta::NRLMSISE00Input, std::shared_ptr< ta::NRLMSISE00Input > >(
            m,
            "NRLMSISE00Input",
            R"doc(Input for computation of NRLMSISE00 atmospheric
                          conditions at current time and position.

                          Input for computation of NRLMSISE00 atmospheric
                          conditions at current time and position. The
                          computation of class may be reperformed every time
                          step, to reflect the changes in atmospheric
                          condition.

                          :param year: Current year
                          :param day_of_year: Day in the current year
                          :param seconds_of_day: Number of seconds into the
                          current day. :param local_solar_time: Local solar
                          time at the computation position :param f107: Current
                          daily F10.7 flux for previous day :param f107a: 81
                          day average of F10.7 flux (centered on current
                          day_of_year). :param ap_daily: Current daily magnetic
                          index :param ap_vector: Current magnetic index data
                          vector: \sa ap_array :param switches: List of
                          NRLMSISE-specific flags: \sa nrlmsise_flags )doc" )
            .def( py::init< int,
                            int,
                            double,
                            double,
                            double,
                            double,
                            double,
                            std::vector< double >,
                            std::vector< int > >( ),
                  py::arg( "year" ) = 0,
                  py::arg( "day_of_year" ) = 0,
                  py::arg( "seconds_of_day" ) = 0.0,
                  py::arg( "local_solar_time" ) = 0.0,
                  py::arg( "f107" ) = 0.0,
                  py::arg( "f107a" ) = 0.0,
                  py::arg( "ap_daily" ) = 0.0,
                  py::arg( "ap_vector" ) = std::vector< double >( 7, 0.0 ),
                  py::arg( "switches" ) = std::vector< int >( ) );

    py::class_< ta::NRLMSISE00Atmosphere, std::shared_ptr< ta::NRLMSISE00Atmosphere > >(
            m,
            "NRLMSISE00Atmosphere",
            R"doc(NRLMSISE00 atmosphere model.

                         This class uses the NRLMSISE00 model to compute the atmospheric density and temperature. The GTD7 function is used: Neutral Atmosphere Empirical Model from the surface to the lower exosphere.

                         Currently, the ideal gas law is used to compute the speed of sound and the specific heat ratio is assumed to be constant and equal to 1.4.

                         :param solar_activity_data: Solar activity data for a range of epochs as produced by tudatpy.io.read_solar_activity_data.
                         )doc" )
            .def( py::init< const std::map<
                          double,
                          std::shared_ptr< tio::solar_activity::SolarActivityData > >, const bool, const bool, const bool >( ),
                  py::arg( "solar_activity_data" ),
                  py::arg( "use_ideal_gas_law" ) = true,
                  py::arg( "use_storm_conditions" ) = false,
                  py::arg( "use_anomalous_oxygen" ) = true )
            .def( "get_density",
                  &ta::NRLMSISE00Atmosphere::getDensity,
                  py::arg( "altitude" ),
                  py::arg( "longitude" ),
                  py::arg( "latitude" ),
                  py::arg( "time" ),
                  R"doc(Get local density

                             Returns the local density at the given altitude,
                             longitude, latitude and time.

                             :param altitude: Altitude at which to get the density. [m]
                             :param longitude: Longitude at which to get the density [rad].
                             :param latitude: Latitude at which to get the density [rad].
                             :param time: Time at which density is to be computed [seconds since J2000].
                             :return: Local density. [kg/m^3]
                             )doc" );

    // END OF NRLMSISE00
    py::enum_< tss::AtmosphereDependentVariables >( m, "AtmosphereDependentVariables" )
            .value( "tabulated_density",
                    tss::AtmosphereDependentVariables::density_dependent_atmosphere )
            .value( "tabulated_pressure",
                    tss::AtmosphereDependentVariables::pressure_dependent_atmosphere )
            .value( "tabulated_temperature",
                    tss::AtmosphereDependentVariables::temperature_dependent_atmosphere )
            .value( "tabulated_gas_constant",
                    tss::AtmosphereDependentVariables::gas_constant_dependent_atmosphere )
            .value( "tabulated_specific_heat_ratio",
                    tss::AtmosphereDependentVariables::specific_heat_ratio_dependent_atmosphere )
            .value( "tabulated_molar_mass",
                    tss::AtmosphereDependentVariables::molar_mass_dependent_atmosphere )
            .export_values( );

    /////////////////////////////////////////////////////////////////////////////
    py::class_< tss::WindModelSettings, std::shared_ptr< tss::WindModelSettings > >(
            m,
            "WindModelSettings",
            R"doc(

         Class for providing settings for wind model.

         Functional (base) class for settings of wind models that require no information in addition to their type.
         Wind model classes requiring additional information must be created using an object derived from this class.





      )doc" );

    py::class_< tss::ConstantWindModelSettings,
                std::shared_ptr< tss::ConstantWindModelSettings >,
                tss::WindModelSettings >(
            m, "ConstantWindModelSettings", R"doc(No documentation found.)doc" );

    py::class_< tss::CustomWindModelSettings,
                std::shared_ptr< tss::CustomWindModelSettings >,
                tss::WindModelSettings >(
            m, "CustomWindModelSettings", R"doc(No documentation found.)doc" );

    py::class_< tss::AtmosphereSettings, std::shared_ptr< tss::AtmosphereSettings > >(
            m,
            "AtmosphereSettings",
            R"doc(

         Base class for providing settings for atmosphere model.

         Functional (base) class for settings of atmosphere models that require no information in addition to their type.
         Atmosphere model classes requiring additional information must be created using an object derived from this class.





      )doc" )
            .def_property( "wind_settings",
                           &tss::AtmosphereSettings::getWindSettings,
                           &tss::AtmosphereSettings::setWindSettings,
                           R"doc(

         **read-only**

         Wind model settings for the atmosphere model settings object.

         :type: WindModelSettings
      )doc" );

    py::class_< tss::ExponentialAtmosphereSettings,
                std::shared_ptr< tss::ExponentialAtmosphereSettings >,
                tss::AtmosphereSettings >( m,
                                           "ExponentialAtmosphereSettings",
                                           R"doc(

         Class for providing settings for exponential atmosphere model.

         `AtmosphereSettings` derived class for a defining the settings of an exponential atmosphere model.




      )doc" );

    py::class_< tss::CustomConstantTemperatureAtmosphereSettings,
                std::shared_ptr< tss::CustomConstantTemperatureAtmosphereSettings >,
                tss::AtmosphereSettings >(
            m, "CustomConstantTemperatureAtmosphereSettings", R"doc(No documentation found.)doc" );

    py::class_< tss::ScaledAtmosphereSettings,
                std::shared_ptr< tss::ScaledAtmosphereSettings >,
                tss::AtmosphereSettings >(
            m, "ScaledAtmosphereSettings", R"doc(No documentation found.)doc" );

    // unexposed this class, because there is no factory
    // function interface yet
    // py::class_<tss::TabulatedAtmosphereSettings,
    //         std::shared_ptr<tss::TabulatedAtmosphereSettings>,
    //         tss::AtmosphereSettings>(m,
    //         "TabulatedAtmosphereSettings",
    //                                  get_docstring("TabulatedAtmosphereSettings").c_str());

    m.def( "constant_wind_model",
           &tss::constantWindModelSettings,
           py::arg( "wind_velocity" ),
           py::arg( "associated_reference_frame" ) = trf::vertical_frame,
           R"doc(

 Function for creating wind model settings with constant wind velocity.

 Function for settings object, defining wind model entirely from constant wind velocity in a given reference frame.


 Parameters
 ----------
 wind_velocity : numpy.ndarray[numpy.float64[3, 1]]
     Constant wind velocity in the specified reference frame.

 associated_reference_frame : numerical_simulation.environment.AerodynamicsReferenceFrames, default = AerodynamicsReferenceFrames.vertical_frame
     Reference frame in which constant wind velocity is defined.

 Returns
 -------
 ConstantWindModelSettings
     Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.WindModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.ConstantWindModelSettings` class





 Examples
 --------
 In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.WindModelSettings`,
 using a constant wind-velocity vector defined in a vertical aerodynamic reference frame:

 .. code-block:: python

   # Define the wind in 3 directions in the vertical reference frame
   wind_Xv = 3     # Meridional wind of +3 m/s (pointing to the North)
   wind_Yv = 5     # Zonal wind of +5 m/s (pointing to the West)
   wind_Zv = -11   # Vertical wind of +11 m/s (pointing out of the centre of the Earth)
   # Create the constant wind settings
   constant_wind = environment_setup.atmosphere.constant_wind_model(
     [wind_Xv, wind_Yv, wind_Zv],
     environment.AerodynamicsReferenceFrames.vertical_frame)
   # Apply the constant wind settings to the Earth atmosphere settings
   body_settings.get("Earth").atmosphere_settings.wind_settings = constant_wind


     )doc" );

    m.def( "custom_wind_model",
           &tss::customWindModelSettings,
           py::arg( "wind_function" ),
           py::arg( "associated_reference_frame" ) = trf::vertical_frame,
           R"doc(

 Function for creating wind model settings with custom wind velocity.

 Function for settings object, defining wind model entirely from custom wind velocity function in a given reference frame.
 The custom wind velocity has to be given as a function of altitude, longitude, latitude and time.

 .. note:: The longitude and latitude will be passed to the function in **degree** and not in radians.
           The altitude is in meters, and the time is a Julian date in seconds since J2000.


 Parameters
 ----------
 wind_velocity : callable[[float, float, float, float], numpy.ndarray[numpy.float64[3, 1]]]
     Custom wind velocity function (w.r.t. altitude, longitude, latitude and time) in the specified reference frame.

 associated_reference_frame : numerical_simulation.environment.AerodynamicsReferenceFrames, default = AerodynamicsReferenceFrames.vertical_frame
     Reference frame in which wind velocity is defined.

 Returns
 -------
 CustomWindModelSettings
     Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.WindModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.CustomWindModelSettings` class





 Examples
 --------
 In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.WindModelSettings`,
 using a user-defined wind-velocity function (of altitude, longitude, latitude and time), defined in a vertical aerodynamic reference frame:

 .. code-block:: python

   # Define the wind in 3 directions in the vertical reference frame
   def wind_function(h, lon, lat, time):
       # Meridional wind (pointing North) depends on latitude [deg] and time [sec since J2000]
       wind_Xv = lat*10/time
       # Zonal wind (pointing West) only depends on the longitude [deg]
       wind_Yv = 5/lon
       # Vertical wind (pointing out of the centre of the Earth) only depends on the altitude [m]
       wind_Zv = 1000/h
       # Return the custom wind
       return [wind_Xv, wind_Yv, wind_Zv]
   # Create the custom wind settings
   custom_wind = environment_setup.atmosphere.custom_wind_model(
       wind_function,
       environment.AerodynamicsReferenceFrames.vertical_frame)
   # Apply the custom wind settings to the Earth atmosphere settings
   body_settings.get("Earth").atmosphere_settings.wind_settings = custom_wind


     )doc" );

    m.def( "exponential_predefined",
           py::overload_cast< const std::string & >( &tss::exponentialAtmosphereSettings ),
           py::arg( "body_name" ),
           R"doc(

 Function for creating atmospheric model settings from pre-defined exponential model.

 Function for settings object, defining atmosphere model from pre-defined exponential model.
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


 Parameters
 ----------
 body_name : str
     Body for which pre-defined model settings are to be loaded. Available bodies "Earth", "Mars".

 Returns
 -------
 ExponentialAtmosphereSettings
     Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.ExponentialAtmosphereSettings` class





 Examples
 --------
 In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` for Mars,
 using the interface of the predefined exponential model, using pre-encoded values:

 .. code-block:: python

    # Create atmosphere settings and add to body settings of "Mars"
    body_settings.get("Mars").atmosphere_settings = environment_setup.atmosphere.exponential_predefined("Mars")


     )doc" );

    m.def( "exponential",
           py::overload_cast< const double,
                              const double,
                              const double,
                              const double,
                              const double >( &tss::exponentialAtmosphereSettings ),
           py::arg( "scale_height" ),
           py::arg( "surface_density" ),
           py::arg( "constant_temperature" ) = 288.15,
           py::arg( "specific_gas_constant" ) =
                   tudat::physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
           py::arg( "ratio_specific_heats" ) = 1.4,
           R"doc(

 Function for creating atmospheric model settings from fully parametrized exponential model.

 Function for settings object, defining exponential atmosphere model.
 The model is solely based on an exponentially decaying density profile with a constant temperature and composition
 (i.e. independent of time, latitude and longitude).

 The user has access to a fully parametrized model, meaning that in addition to the required input parameters ``scale_height`` and ``surface_density`` (ground-level air density),
 the user can specify non-standard values for constant temperature, gas constant and specific heats ratio.


 Parameters
 ----------
 scale_height : float
     Scale height for density profile of atmosphere.
 surface_density : float
     Atmospheric density at ground level.
 constant_temperature : float, default = 288.15
     Constant atmospheric temperature.
 specific_gas_constant : float, default = constants.SPECIFIC_GAS_CONSTANT_AIR
     Specific gas constant for (constant) atmospheric chemical composition.
 ratio_specific_heats : float, default = 1.4
     Ratio of specific heats for (constant) atmospheric chemical composition.
 Returns
 -------
 ExponentialAtmosphereSettings
     Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.ExponentialAtmosphereSettings` class





 Examples
 --------
 In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` for Earth,
 using the minimalist interface to the exponential model and taking parameters with classic values for Earth:

 .. code-block:: python

    # define parameters of an invariant exponential atmosphere model
    density_scale_height = 7.2E3
    constant_temperature = 290
    # create atmosphere settings and add to body settings of "Earth"
    body_settings.get( "Earth" ).atmosphere_settings = environment_setup.atmosphere.exponential(
        density_scale_height, density_at_zero_altitude)


     )doc" );

    m.def( "nrlmsise00",
           &tss::nrlmsise00AtmosphereSettings,
           py::arg( "space_weather_file" ) =
                   tudat::paths::getSpaceWeatherDataPath( ) + "/sw19571001.txt",
           py::arg( "use_storm_conditions" ) = false,
           py::arg( "use_anomalous_oxygen" ) = true,
               R"doc(

 Function for creating NRLMSISE-00 atmospheric model settings.

 Function for settings object, defining atmosphere model in accordance to the NRLMSISE-00 global reference model for Earth's atmosphere.
 The NRLMSISE-00 model implementation uses the code `here <https://github.com/tudat-team/nrlmsise-00-cmake>`_.


 Parameters
 ----------
space_weather_file : str, default = :func:`~tudatpy.data.get_space_weather_path` + 'sw19571001.txt'
    File to be used for space weather characteristics as a function of time (e.g. F10.7, Kp, etc.). The file is typically taken from here `celestrak <https://celestrak.org/SpaceData/sw19571001.txt>`_ (note that the file in your resources path will not be the latest version of this file; download and replace your existing file if required). Documentation on the file is given `here <https://celestrak.org/SpaceData/SpaceWx-format.php>`_
use_storm_conditions : bool, default = false
    Boolean to define whether to use sub-daily Ap values when querying the NRLMSISE model, which is relevant under geomagnetic storm conditions (see `NRLMSISE code <https://github.com/tudat-team/nrlmsise-00-cmake/blob/master/nrlmsise-00.h>`_, setting this variable to true sets ``switches[9]`` to -1, with resulting details of Ap values defined in ``ap_array``).
use_anomalous_oxygen : bool, default = true
    Boolean to define whether to use anomalous oxygen when querying the NRLMSISE model (if true, using ``gtd7d`` function, if false using ``gtd7`` function in `NRLMSISE code <https://github.com/tudat-team/nrlmsise-00-cmake/blob/master/nrlmsise-00.h>`_)
Returns
 -------
 AtmosphereSettings
     Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` class





 Examples
 --------
 In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` for Earth,
 using the NRLMSISE-00 global reference model:

 .. code-block:: python

    # create atmosphere settings and add to body settings of body "Earth"
    body_settings.get( "Earth" ).atmosphere_settings = environment_setup.atmosphere.nrlmsise00()


     )doc" );

    m.def( "tabulated",
           &tss::tabulatedAtmosphereSettings,
           py::arg( "atmosphere_data_file" ),
           py::arg( "dependent_variable_names" ) = std::vector< ta::AtmosphereDependentVariables >(
                   { tss::density_dependent_atmosphere,
                     tss::pressure_dependent_atmosphere,
                     tss::temperature_dependent_atmosphere } ),
           py::arg( "specific_gas_constant" ) = tp::SPECIFIC_GAS_CONSTANT_AIR,
           py::arg( "ratio_of_specific_heats" ) = 1.4 );

    m.def( "us76",
           &tss::us76AtmosphereSettings,
           R"doc(

 Function for creating US76 standard atmosphere model settings.

 Function for creating US76 standard atmosphere model settings. The model is defined using tabulated data for density, pressure and temperature,
 from an altitude of -5 km up to 1000 km. Up to 100 km, a data point is provided every 100 m. Above 100 km, a data point is provided every 1 km. The data
 are interpolated using a cubic spline interpolator. Note that this model is specific to Earth's atmosphere.

 Returns
 -------
 AtmosphereSettings
     Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` class





 Examples
 --------
 In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` for Earth,
 using the US76 standard atmosphere model:

 .. code-block:: python

    # create atmosphere settings and add to body settings of body "Earth"
    body_settings.get( "Earth" ).atmosphere_settings = environment_setup.atmosphere.us76()




     )doc" );

    m.def( "custom_constant_temperature",
           py::overload_cast< const std::function< double( const double ) >,
                              const double,
                              const double,
                              const double >( &tss::customConstantTemperatureAtmosphereSettings ),
           py::arg( "density_function" ),
           py::arg( "constant_temperature" ),
           py::arg( "specific_gas_constant" ) =
                   tudat::physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
           py::arg( "ratio_of_specific_heats" ) = 1.4,
           R"doc(

 Function for creating atmospheric model settings from custom density profile.

 Function for settings object, defining constant temperature atmosphere model from custom density profile.
 The user is specifying the density profile as a function of altitude.
 The value of pressure is computed by assuming hydrostatic equilibrium, temperature, gas constant and the ratio of specific heats are modelled as constants.


 Parameters
 ----------
 density_function : callable[[float], float]
     Function to retrieve the density at the current altitude.

 constant_temperature : float
     Constant atmospheric temperature.
 specific_gas_constant : float, default = 287.0
     Specific gas constant for (constant) atmospheric chemical composition.
 ratio_specific_heats : float, default = 1.4
     Ratio of specific heats for (constant) atmospheric chemical composition.
 Returns
 -------
 CustomConstantTemperatureAtmosphereSettings
     Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.CustomConstantTemperatureAtmosphereSettings` class





 Examples
 --------
 In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` for Earth,
 with constant temperature and composition, but a density which varies with altitude according to a user-defined model:

 .. code-block:: python

   # Define the density as a function of altitude [in m]
   def density_function(h):
       # Return the density according to a modified exponential model
       return 1.15 * np.exp(-h/7300)
   # Define parameters for constant temperature and composition
   constant_temperature = 250.0
   specific_gas_constant = 300.0
   ratio_of_specific_heats = 1.4
   # Create the custom constant temperature atmosphere settings
   custom_density_settings = environment_setup.atmosphere.custom_constant_temperature(
       density_function,
       constant_temperature,
       specific_gas_constant,
       ratio_of_specific_heats)
   # Add the custom density to the body settings of "Earth"
   body_settings.get("Earth").atmosphere_settings = custom_density_settings


     )doc" );

    m.def( "custom_four_dimensional_constant_temperature",
           py::overload_cast< const std::function< double(
                                      const double, const double, const double, const double ) >,
                              const double,
                              const double,
                              const double >( &tss::customConstantTemperatureAtmosphereSettings ),
           py::arg( "density_function" ),
           py::arg( "constant_temperature" ),
           py::arg( "specific_gas_constant" ) =
                   tudat::physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
           py::arg( "ratio_of_specific_heats" ) = 1.4,
           R"doc(

 Function for creating atmospheric model settings from custom density profile.

 Function for settings object, defining constant temperature atmosphere model from custom density profile.
 The user is specifying the density profile as a function of altitude, longitude, latitude and time.

 .. note:: The longitude and latitude will be passed to the function in **degree** and not in radians.
           The altitude is in meters, and the time is a Julian date in seconds since J2000.


 Parameters
 ----------
 density_function : callable[[float, float, float, float], float]
     Function to retrieve the density at the current altitude, longitude, latitude and time.

 constant_temperature : float
     Constant atmospheric temperature.
 specific_gas_constant : float, default = 287.0
     Specific gas constant for (constant) atmospheric chemical composition.
 ratio_specific_heats : float, default = 1.4
     Ratio of specific heats for (constant) atmospheric chemical composition.
 Returns
 -------
 CustomConstantTemperatureAtmosphereSettings
     Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.CustomConstantTemperatureAtmosphereSettings` class





 Examples
 --------
 In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` for Earth,
 with constant temperature and composition (gas constant and ratio of specific heats), but a density which varies with altitude, longitude, latitude and time, according to a user-defined model:

 .. code-block:: python

   # Define the density as a function of altitude [m], longitude [deg], latitude [deg], and time [sec since J2000]
   def density_function(h, lon, lat, time):
       # Return the density according to an exponential model that varies with time to add noise with a sine (ignore lon/lat)
       return (1 + 0.15 * np.sin(time/10)) * np.exp(-h/7300)
   # Define the parameters for constant temperature and composition
   constant_temperature = 250.0
   specific_gas_constant = 300.0
   ratio_of_specific_heats = 1.4
   # Create the atmosphere settings and add to body settings of "Earth"
   body_settings.get( "Earth" ).atmosphere_settings = environment_setup.atmosphere.custom_constant_temperature(
       density_function,
       constant_temperature,
       specific_gas_constant,
       ratio_of_specific_heats )


     )doc" );

    m.def( "scaled_by_function",
           py::overload_cast< const std::shared_ptr< tss::AtmosphereSettings >,
                              const std::function< double( const double ) >,
                              const bool >( &tss::scaledAtmosphereSettings ),
           py::arg( "unscaled_atmosphere_settings" ),
           py::arg( "density_scaling_function" ),
           py::arg( "is_scaling_absolute" ) = false,
           R"doc(

 Function for creating scaled atmospheric model settings.

 Function for settings object, defining atmospheric model based on scaling an existing atmospheric settings object.
 The user can apply custom scaling factors (or absolute values) to the air densities of the existing model settings (for instance for an uncertainty analysis).


 Parameters
 ----------
 unscaled_atmosphere_settings : AtmosphereSettings
     Sets base settings of atmosphere model to be scaled.
 density_scaling_function : Callable[[float], float]
     Specifies air density scaling factor as a function of time.
 is_scaling_absolute : bool, default=false
     Boolean indicating whether density scaling is absolute. Setting this boolean to true will add the scaling value to the baseline density, instead of the default behaviour of multiplying the baseline density by the scaling value.

 Returns
 -------
 ScaledAtmosphereSettings
     Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.ScaledAtmosphereSettings` class.



 Notes
 -----
 At present, the scaled atmosphere model only supports scaling of the density value.
 For cases where the density is used to compute other atmospheric quantities (such as pressure using hydrostatic equilibrium),
 this calculation is performed using the `unscaled` density!



 Examples
 --------
 In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` for Earth,
 by modifying an existing :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` object such, that the resulting air density profile is scaled with a user-defined function of time:

 .. code-block:: python

   # Define the density scaling as a function of time [sec since J2000] (to add noise with a sine)
   def scaling_function(time):
       return 1 + np.sin(time / 50) * 0.25
   # Extract the existing atmosphere model settings
   unscaled_atmosphere_settings = body_settings.get( "Earth" ).atmosphere_settings
   # Create the atmosphere settings and add to body settings of "Earth"
   body_settings.get( "Earth" ).atmosphere_settings =  environment_setup.atmosphere.scaled_by_function(
       unscaled_atmosphere_settings,
       scaling_function )


     )doc" );

    m.def( "scaled_by_constant",
           py::overload_cast< const std::shared_ptr< tss::AtmosphereSettings >,
                              const double,
                              const bool >( &tss::scaledAtmosphereSettings ),
           py::arg( "unscaled_atmosphere_settings" ),
           py::arg( "density_scaling" ),
           py::arg( "is_scaling_absolute" ) = false,
           R"doc(

 Function for creating scaled atmospheric model settings.

 Function for settings object, defining atmospheric model based on an scaling of an existing atmospheric settings object.
 The user can apply a scaling factor (or an absolute value) to the air densities of the existing model settings (for instance for an uncertainty analysis).


 Parameters
 ----------
 unscaled_atmosphere_settings : AtmosphereSettings
     Sets base settings of atmosphere model to be scaled.
 density_scaling : float
     Constant scaling factor to be applied to the entire air density profile.
 is_scaling_absolute : bool, default=false
     Boolean indicating whether density scaling is absolute. Setting this boolean to true will add the scaling value to the baseline density, instead of the default behaviour of multiplying the baseline density by the scaling value.

 Returns
 -------
 ScaledAtmosphereSettings
     Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.ScaledAtmosphereSettings` class.



 Notes
 -----
 At present, the scaled atmosphere model only supports scaling of the density value.
 For cases where the density is used to compute other atmospheric quantities (such as pressure using hydrostatic equilibrium),
 this calculation is performed using the `unscaled` density!



 Examples
 --------
 In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` for Earth,
 by modifying an existing :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` object such that the resulting air density profile is scaled by a constant:

 .. code-block:: python

   # define parameter for scaling
   scaling_constant = 1.5
   # define variable containing the existing atmosphere model settings
   unscaled_atmosphere_settings = body_settings.get( "Earth" ).atmosphere_settings
   # create atmosphere settings and add to body settings of "Earth"
   body_settings.get( "Earth" ).atmosphere_settings =  environment_setup.atmosphere.scaled_by_constant(
       unscaled_atmosphere_settings,
       scaling_constant )


     )doc" );
}

}  // namespace atmosphere
}  // namespace environment_setup
}  // namespace numerical_simulation
}  // namespace tudatpy
