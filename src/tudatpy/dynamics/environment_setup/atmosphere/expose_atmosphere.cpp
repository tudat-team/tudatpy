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
namespace dynamics
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
            .def( "set_use_geodetic_latitude", &ta::NRLMSISE00Atmosphere::setUseGeodeticLatitude)
            .def( "get_use_geodetic_latitude", &ta::NRLMSISE00Atmosphere::getUseGeodeticLatitude)
            .def( "set_use_utc", &ta::NRLMSISE00Atmosphere::setUseUtc)
            .def( "get_use_utc", &ta::NRLMSISE00Atmosphere::getUseUtc)
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





      )doc" )
            .def_property( "include_corotation",
                           &tss::WindModelSettings::getIncludeCorotation,
                           &tss::WindModelSettings::setIncludeCorotation,
                           R"doc(

         Boolean flag indicating whether atmospheric co-rotation should be included in aerodynamic computations.

         :type: bool
      )doc" );

    py::class_< tss::EmptyWindModelSettings,
                std::shared_ptr< tss::EmptyWindModelSettings >,
                tss::WindModelSettings >(
            m, "EmptyWindModelSettings", R"doc(Settings for empty wind model (no physical wind, only co-rotation control).)doc" );

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

    m.def( "empty_wind_model",
           &tss::emptyWindModelSettings,
           py::arg( "include_corotation" ) = true,
           R"doc(

 Function for creating empty wind model settings.

 Function for settings object for an empty wind model (no physical wind, returns zero velocity).
 This is useful when you want to control atmospheric co-rotation behavior without specifying actual wind.


 Parameters
 ----------
 include_corotation : bool, default = True
     Boolean flag indicating whether atmospheric co-rotation should be included in aerodynamic computations.

 Returns
 -------
 EmptyWindModelSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.WindModelSettings` derived :class:`~tudatpy.dynamics.environment_setup.atmosphere.EmptyWindModelSettings` class


 Examples
 --------
 In this example, we create :class:`~tudatpy.dynamics.environment_setup.atmosphere.WindModelSettings`,
 for an atmosphere without physical wind but with co-rotation disabled:

 .. code-block:: python

   # Create empty wind model with co-rotation disabled
   empty_wind = environment_setup.atmosphere.empty_wind_model(include_corotation=False)
   # Apply to the atmosphere settings
   body_settings.get("Earth").atmosphere_settings.wind_settings = empty_wind


     )doc" );

    m.def( "constant_wind_model",
           &tss::constantWindModelSettings,
           py::arg( "wind_velocity" ),
           py::arg( "associated_reference_frame" ) = trf::vertical_frame,
           py::arg( "include_corotation" ) = true,
           R"doc(

 Function for creating wind model settings with constant wind velocity.

 Function for settings object, defining wind model entirely from constant wind velocity in a given reference frame.


 Parameters
 ----------
 wind_velocity : numpy.ndarray[numpy.float64[3, 1]]
     Constant wind velocity in the specified reference frame.

 associated_reference_frame : dynamics.environment.AerodynamicsReferenceFrames, default = AerodynamicsReferenceFrames.vertical_frame
     Reference frame in which constant wind velocity is defined.

 include_corotation : bool, default = True
     Boolean flag indicating whether atmospheric co-rotation should be included in aerodynamic computations.

 Returns
 -------
 ConstantWindModelSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.WindModelSettings` derived :class:`~tudatpy.dynamics.environment_setup.atmosphere.ConstantWindModelSettings` class





 Examples
 --------
 In this example, we create :class:`~tudatpy.dynamics.environment_setup.atmosphere.WindModelSettings`,
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
           py::arg( "include_corotation" ) = true,
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

 associated_reference_frame : dynamics.environment.AerodynamicsReferenceFrames, default = AerodynamicsReferenceFrames.vertical_frame
     Reference frame in which wind velocity is defined.

 include_corotation : bool, default = True
     Boolean flag indicating whether atmospheric co-rotation should be included in aerodynamic computations.

 Returns
 -------
 CustomWindModelSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.WindModelSettings` derived :class:`~tudatpy.dynamics.environment_setup.atmosphere.CustomWindModelSettings` class





 Examples
 --------
 In this example, we create :class:`~tudatpy.dynamics.environment_setup.atmosphere.WindModelSettings`,
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

    m.def( "coma_wind_model",
           &tss::comaWindModelSettings,
           py::arg( "dataset_collection" ),
           py::arg( "requested_max_degree" ) = -1,
           py::arg( "requested_max_order" ) = -1,
           py::arg( "associated_reference_frame" ) = trf::vertical_frame,
           py::arg( "include_corotation" ) = true,
           R"doc(

 Function for creating coma wind model settings.

 Function for settings object, defining coma wind model from a dataset collection containing
 x, y, z wind velocity components. The wind model uses spherical harmonic expansion to compute
 wind velocities as a function of position.


 Parameters
 ----------
 dataset_collection : ComaWindDatasetCollection
     Collection containing x, y, z component datasets (either polynomial or Stokes coefficients).

 requested_max_degree : int, default = -1
     Maximum spherical harmonic degree to use (-1 for automatic determination from data).

 requested_max_order : int, default = -1
     Maximum spherical harmonic order to use (-1 for automatic determination from data).

 associated_reference_frame : dynamics.environment.AerodynamicsReferenceFrames, default = AerodynamicsReferenceFrames.vertical_frame
     Reference frame in which the wind velocity is defined.

 include_corotation : bool, default = True
     Boolean flag indicating whether atmospheric co-rotation should be included in aerodynamic computations.

 Returns
 -------
 WindModelSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.WindModelSettings` derived :class:`~tudatpy.dynamics.environment_setup.atmosphere.ComaWindModelSettings` class


 Examples
 --------
 In this example, we create :class:`~tudatpy.dynamics.environment_setup.atmosphere.WindModelSettings`
 for a coma wind model using a dataset collection:

 .. code-block:: python

   # Create file processor from polynomial coefficient files
   wind_processor = environment_setup.atmosphere.coma_wind_file_processor(
       x_file_paths, y_file_paths, z_file_paths)

   # Create dataset collection with Stokes coefficients
   wind_datasets = wind_processor.create_sh_datasets(
       radii_m=[1000, 2000, 3000],
       sol_longitudes_deg=[0, 90, 180, 270])

   # Create coma wind model settings in vertical frame
   coma_wind = environment_setup.atmosphere.coma_wind_model(
       wind_datasets,
       requested_max_degree=10,
       requested_max_order=10,
       associated_reference_frame=environment.AerodynamicsReferenceFrames.vertical_frame,
       include_corotation=True)

   # Apply to atmosphere settings
   body_settings.get("Comet").atmosphere_settings.wind_settings = coma_wind


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
     Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.dynamics.environment_setup.atmosphere.ExponentialAtmosphereSettings` class





 Examples
 --------
 In this example, we create :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` for Mars,
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
     Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.dynamics.environment_setup.atmosphere.ExponentialAtmosphereSettings` class





 Examples
 --------
 In this example, we create :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` for Earth,
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
The NRLMSISE-00 model implementation uses the code from `tudat-team/nrlmsise-00-cmake <https://github.com/tudat-team/nrlmsise-00-cmake>`_.


Parameters
----------
space_weather_file : str, default = :func:`~tudatpy.data.get_space_weather_path` + 'sw19571001.txt'
    File to be used for space weather characteristics as a function of time (e.g. F10.7, Kp, etc.). The file is typically taken from `celestrak <https://celestrak.org/SpaceData/sw19571001.txt>`_ (note that the file in your resources path will not be the latest version of this file; download and replace your existing file if required). Documentation on the file is given on the `celestrak website <https://celestrak.org/SpaceData/SpaceWx-format.php>`_
use_storm_conditions : bool, default = false
    Boolean to define whether to use sub-daily Ap values when querying the NRLMSISE model, which is relevant under geomagnetic storm conditions (see `NRLMSISE code <https://github.com/tudat-team/nrlmsise-00-cmake/blob/master/nrlmsise-00.h>`_, setting this variable to true sets ``switches[9]`` to -1, with resulting details of Ap values defined in ``ap_array``).
use_anomalous_oxygen : bool, default = true
    Boolean to define whether to use anomalous oxygen when querying the NRLMSISE model (if true, using ``gtd7d`` function, if false using ``gtd7`` function in `NRLMSISE code <https://github.com/tudat-team/nrlmsise-00-cmake/blob/master/nrlmsise-00.h>`_)

Returns
-------
AtmosphereSettings
    Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` for Earth,
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
     Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` class





 Examples
 --------
 In this example, we create :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` for Earth,
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
     Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.dynamics.environment_setup.atmosphere.CustomConstantTemperatureAtmosphereSettings` class





 Examples
 --------
 In this example, we create :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` for Earth,
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
     Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.dynamics.environment_setup.atmosphere.CustomConstantTemperatureAtmosphereSettings` class





 Examples
 --------
 In this example, we create :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` for Earth,
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
     Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.dynamics.environment_setup.atmosphere.ScaledAtmosphereSettings` class.



 Notes
 -----
 At present, the scaled atmosphere model only supports scaling of the density value.
 For cases where the density is used to compute other atmospheric quantities (such as pressure using hydrostatic equilibrium),
 this calculation is performed using the `unscaled` density!



 Examples
 --------
 In this example, we create :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` for Earth,
 by modifying an existing :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` object such, that the resulting air density profile is scaled with a user-defined function of time:

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
     Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.dynamics.environment_setup.atmosphere.ScaledAtmosphereSettings` class.



 Notes
 -----
 At present, the scaled atmosphere model only supports scaling of the density value.
 For cases where the density is used to compute other atmospheric quantities (such as pressure using hydrostatic equilibrium),
 this calculation is performed using the `unscaled` density!



 Examples
 --------
 In this example, we create :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` for Earth,
 by modifying an existing :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` object such that the resulting air density profile is scaled by a constant:

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

    // --- Coma Model ---

    m.def(
            "coma_model",
            py::overload_cast<
                const tss::ComaPolyDataset&, double, int, int >( &tss::comaSettings ),
            py::arg( "poly_data" ),
            py::arg( "molecular_weight" ),
            py::arg( "max_degree" ) = -1,
            py::arg( "max_order" ) = -1,
            R"doc(

 Function for creating coma atmosphere model settings from polynomial coefficients.

 Function for settings object, defining a coma atmosphere model based on spherical harmonic expansion
 of gas density data. The coma model is designed for modeling cometary atmospheres (comae) where
 gas density varies with position and time. This variant uses polynomial coefficient data that
 describes the spatial distribution of gas density.

 The density is computed using spherical harmonic expansion, allowing efficient representation of
 complex 3D density distributions around the nucleus. The model supports time-dependent density
 variations through multiple data files covering different time periods.


 Parameters
 ----------
 poly_data : ComaPolyDataset
     Polynomial coefficient dataset containing spherical harmonic coefficients for gas density
     distribution. Create using :func:`~tudatpy.dynamics.environment_setup.atmosphere.coma_model_file_processor`.

 molecular_weight : float
     Molecular weight (molar mass) of the gas species [kg/mol]. For water vapor (H2O), use 0.018015 kg/mol.

 max_degree : int, default = -1
     Maximum spherical harmonic degree to use in density calculations. Set to -1 to automatically
     use the maximum degree available in the dataset.

 max_order : int, default = -1
     Maximum spherical harmonic order to use in density calculations. Set to -1 to automatically
     use the maximum order available in the dataset.

 Returns
 -------
 AtmosphereSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` class
     configured for coma model.


 Examples
 --------
 In this example, we create a coma atmosphere model from polynomial coefficient files:

 .. code-block:: python

   # Define paths to polynomial coefficient files
   poly_file_paths = [
       "coma_data/poly_coeffs_epoch1.txt",
       "coma_data/poly_coeffs_epoch2.txt"
   ]

   # Create file processor from polynomial files
   processor = environment_setup.atmosphere.coma_model_file_processor(poly_file_paths)

   # Create polynomial dataset
   poly_dataset = processor.create_poly_coef_dataset()

   # Create coma atmosphere settings
   coma_atmosphere = environment_setup.atmosphere.coma_model(
       poly_data=poly_dataset,
       molecular_weight=0.018015,  # H2O molecular weight in kg/mol
       max_degree=10,
       max_order=10)

   # Apply to body settings
   body_settings.get("67P").atmosphere_settings = coma_atmosphere


    )doc"
                );

    m.def(
            "coma_model",
            py::overload_cast<
                const tss::ComaStokesDataset&, double, int, int >( &tss::comaSettings ),
            py::arg( "stokes_data" ),
            py::arg( "molecular_weight" ),
            py::arg( "max_degree" ) = -1,
            py::arg( "max_order" ) = -1,
            R"doc(

 Function for creating coma atmosphere model settings from Stokes coefficients.

 Function for settings object, defining a coma atmosphere model based on spherical harmonic expansion
 of gas density data using precomputed Stokes coefficients. The coma model is designed for modeling
 cometary atmospheres (comae) where gas density varies with position and time. This variant uses
 precomputed Stokes coefficients (spherical harmonics) evaluated at specific radii and solar longitudes.

 Stokes coefficients provide a more direct representation of the spherical harmonic expansion compared
 to polynomial coefficients, offering faster evaluation during simulation. The coefficients
 are pre-evaluated at a grid of radii and solar longitudes, with interpolation used for intermediate values.


 Parameters
 ----------
 stokes_data : ComaStokesDataset
     Precomputed Stokes coefficient dataset containing spherical harmonic coefficients evaluated at
     specific radii and solar longitudes. Create using :func:`~tudatpy.dynamics.environment_setup.atmosphere.coma_model_file_processor`
     or load from pre-existing Stokes coefficient CSV files.

 molecular_weight : float
     Molecular weight (molar mass) of the gas species [kg/mol]. For water vapor (H2O), use 0.018015 kg/mol.

 max_degree : int, default = -1
     Maximum spherical harmonic degree to use in density calculations. Set to -1 to automatically
     use the maximum degree available in the dataset.

 max_order : int, default = -1
     Maximum spherical harmonic order to use in density calculations. Set to -1 to automatically
     use the maximum order available in the dataset.

 Returns
 -------
 AtmosphereSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` class
     configured for coma model.


 Examples
 --------
 In this example, we create a coma atmosphere model by converting polynomial coefficients to Stokes coefficients:

 .. code-block:: python

   # Create file processor from polynomial coefficient files
   poly_file_paths = ["coma_data/poly_coeffs_epoch1.txt"]
   processor = environment_setup.atmosphere.coma_model_file_processor(poly_file_paths)

   # Create Stokes dataset by evaluating at specific radii and solar longitudes
   stokes_dataset = processor.create_sh_dataset(
       radii_m=[1000.0, 2000.0, 5000.0, 10000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0],
       requested_max_degree=10,
       requested_max_order=10)

   # Create coma atmosphere settings
   coma_atmosphere = environment_setup.atmosphere.coma_model(
       stokes_data=stokes_dataset,
       molecular_weight=0.018015,  # H2O molecular weight in kg/mol
       max_degree=10,
       max_order=10)

   # Apply to body settings
   body_settings.get("67P").atmosphere_settings = coma_atmosphere

 Alternatively, load from pre-existing Stokes coefficient CSV files:

 .. code-block:: python

   # Create file processor from existing Stokes CSV files
   processor = environment_setup.atmosphere.coma_model_file_processor(
       input_dir="coma_data/stokes_files",
       prefix="stokes")

   # Create Stokes dataset (radii and longitudes are read from files)
   stokes_dataset = processor.create_sh_dataset(
       radii_m=[],  # Ignored when loading from files
       sol_longitudes_deg=[])

   # Create coma atmosphere settings
   coma_atmosphere = environment_setup.atmosphere.coma_model(
       stokes_data=stokes_dataset,
       molecular_weight=0.018015)

   # Apply to body settings
   body_settings.get("67P").atmosphere_settings = coma_atmosphere


    )doc"
            );

    // === Coma processing: datasets (minimal shells so Python can hold them) ===
    py::class_< tss::ComaPolyDataset >( m,
                                        "ComaPolyDataset",
                                        R"doc(

 Polynomial coefficient dataset for coma atmosphere density model.

 This class holds polynomial coefficients that describe the spherical harmonic expansion of gas
 density in a cometary coma. The coefficients represent the spatial variation of density as a
 function of position relative to the comet nucleus. The dataset can contain data from multiple
 files, each covering a different time period, enabling time-dependent coma modeling.

 Polynomial coefficients are the raw input format and need to be evaluated to produce Stokes
 coefficients for use in the coma model. This evaluation is typically performed automatically
 when creating a coma atmosphere model.

 .. note:: This class cannot be directly instantiated by the user. Create instances using
           :func:`~tudatpy.dynamics.environment_setup.atmosphere.coma_model_file_processor`
           and its ``create_poly_coef_dataset()`` method.


 Examples
 --------
 Create a polynomial dataset from coefficient files:

 .. code-block:: python

   # Define paths to polynomial coefficient files
   file_paths = [
       "coma_data/h2o_poly_epoch1.txt",
       "coma_data/h2o_poly_epoch2.txt"
   ]

   # Create file processor
   processor = environment_setup.atmosphere.coma_model_file_processor(file_paths)

   # Create polynomial dataset
   poly_dataset = processor.create_poly_coef_dataset()

   # Use the dataset to create coma atmosphere settings
   coma_settings = environment_setup.atmosphere.coma_model(
       poly_data=poly_dataset,
       molecular_weight=0.018015)  # H2O molecular weight


      )doc" );

    py::class_< tss::ComaStokesDataset >( m,
                                          "ComaStokesDataset",
                                          R"doc(

 Stokes coefficient (spherical harmonics) dataset for coma atmosphere density model.

 This class holds precomputed Stokes coefficients (spherical harmonic coefficients) that describe
 the gas density distribution in a cometary coma. The coefficients are evaluated at a specific grid
 of radii and solar longitudes, enabling efficient interpolation during simulation. The dataset can
 contain data from multiple files, each covering a different time period.

 Stokes coefficients provide a direct representation of the spherical harmonic expansion:

 .. math::
     \\rho(r, \\theta, \\phi) = \\sum_{n=0}^{N} \\sum_{m=0}^{n} [C_{nm}(r) \\cos(m\\phi) + S_{nm}(r) \\sin(m\\phi)] P_{nm}(\\cos\\theta)

 where :math:`C_{nm}` and :math:`S_{nm}` are the Stokes coefficients, :math:`P_{nm}` are associated
 Legendre polynomials, and :math:`(r, \\theta, \\phi)` are spherical coordinates.

 .. note:: This class cannot be directly instantiated by the user. Create instances using
           :func:`~tudatpy.dynamics.environment_setup.atmosphere.coma_model_file_processor`
           and its ``create_sh_dataset()`` method, or load from pre-existing CSV files.


 Examples
 --------
 Create a Stokes dataset by transforming polynomial coefficients:

 .. code-block:: python

   # Create file processor from polynomial files
   processor = environment_setup.atmosphere.coma_model_file_processor(
       ["coma_data/h2o_poly.txt"])

   # Transform to Stokes coefficients at specific radii and solar longitudes
   stokes_dataset = processor.create_sh_dataset(
       radii_m=[1000.0, 2000.0, 5000.0, 10000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0],
       requested_max_degree=10,
       requested_max_order=10)

   # Use the dataset to create coma atmosphere settings
   coma_settings = environment_setup.atmosphere.coma_model(
       stokes_data=stokes_dataset,
       molecular_weight=0.018015)

 Alternatively, load from pre-existing Stokes CSV files:

 .. code-block:: python

   # Create file processor from existing Stokes CSV files
   processor = environment_setup.atmosphere.coma_model_file_processor(
       input_dir="coma_data/stokes_csv",
       prefix="stokes")

   # Load Stokes dataset from files
   stokes_dataset = processor.create_sh_dataset(
       radii_m=[],  # Ignored when loading from files
       sol_longitudes_deg=[])

   # Use the dataset
   coma_settings = environment_setup.atmosphere.coma_model(
       stokes_data=stokes_dataset,
       molecular_weight=0.018015)


      )doc" );

    py::class_< tss::ComaWindDatasetCollection >( m,
                                                   "ComaWindDatasetCollection",
                                                   R"doc(
        Collection of three coma datasets for wind model (x, y, z components).

        This class holds three datasets (one for each spatial component of the wind velocity)
        that are used together to construct a ComaWindModel. All three datasets must be of
        the same type (either all polynomial or all Stokes coefficients).
        )doc" )
        .def("is_poly", &tss::ComaWindDatasetCollection::isPoly,
             R"doc(

 Check if the collection contains polynomial coefficient datasets.

 Returns
 -------
 bool
     True if the collection contains :class:`~tudatpy.dynamics.environment_setup.atmosphere.ComaPolyDataset` objects,
     False if it contains :class:`~tudatpy.dynamics.environment_setup.atmosphere.ComaStokesDataset` objects.


 Examples
 --------
 .. code-block:: python

   # Create wind dataset collection from polynomial files
   processor = environment_setup.atmosphere.coma_wind_file_processor(
       x_file_paths, y_file_paths, z_file_paths)
   wind_datasets = processor.create_poly_coef_datasets()

   # Check dataset type
   if wind_datasets.is_poly():
       print("Collection contains polynomial datasets")


      )doc")
        .def("is_stokes", &tss::ComaWindDatasetCollection::isStokes,
             R"doc(

 Check if the collection contains Stokes coefficient datasets.

 Returns
 -------
 bool
     True if the collection contains :class:`~tudatpy.dynamics.environment_setup.atmosphere.ComaStokesDataset` objects,
     False if it contains :class:`~tudatpy.dynamics.environment_setup.atmosphere.ComaPolyDataset` objects.


 Examples
 --------
 .. code-block:: python

   # Create wind dataset collection from polynomial files and convert to Stokes
   processor = environment_setup.atmosphere.coma_wind_file_processor(
       x_file_paths, y_file_paths, z_file_paths)
   wind_datasets = processor.create_sh_datasets(
       radii_m=[1000.0, 2000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0])

   # Check dataset type
   if wind_datasets.is_stokes():
       print("Collection contains Stokes datasets")


      )doc");

    // === Coma processing: file processor ===
    py::class_< tss::ComaModelFileProcessor, std::shared_ptr< tss::ComaModelFileProcessor > >( m,
                                               "ComaModelFileProcessor",
                                               R"doc(

 Processor for creating coma atmosphere datasets from coefficient files.

 This class provides a high-level interface for loading and processing coma model data from files
 containing either polynomial coefficients or pre-computed Stokes (spherical harmonic) coefficients.
 It handles two distinct workflows:

 1. **Polynomial coefficient workflow**: Load polynomial coefficients from text files, optionally
    transform them to Stokes coefficients at specified radii and solar longitudes, and create datasets
    for use in coma atmosphere models.

 2. **Stokes coefficient workflow**: Load pre-computed Stokes coefficients from CSV files that were
    previously generated and saved.

 The processor automatically manages file reading, data validation, and coordinate transformations,
 simplifying the setup of complex coma atmosphere models.

 .. note:: This class cannot be directly instantiated. Create instances using the factory functions
           :func:`~tudatpy.dynamics.environment_setup.atmosphere.coma_model_file_processor`.


 Examples
 --------
 Polynomial coefficient workflow:

 .. code-block:: python

   # Create processor from polynomial coefficient files
   poly_files = ["h2o_epoch1.txt", "h2o_epoch2.txt"]
   processor = environment_setup.atmosphere.coma_model_file_processor(poly_files)

   # Create polynomial dataset directly
   poly_dataset = processor.create_poly_coef_dataset()

   # Or transform to Stokes coefficients
   stokes_dataset = processor.create_sh_dataset(
       radii_m=[1000.0, 2000.0, 5000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0])

 Stokes coefficient workflow:

 .. code-block:: python

   # Create processor from existing Stokes CSV files
   processor = environment_setup.atmosphere.coma_model_file_processor(
       input_dir="coma_stokes_data",
       prefix="stokes")

   # Load Stokes dataset
   stokes_dataset = processor.create_sh_dataset(
       radii_m=[],  # Ignored when loading from files
       sol_longitudes_deg=[])


      )doc" )
        .def("create_poly_coef_dataset",
             &tss::ComaModelFileProcessor::createPolyCoefDataset,
             R"doc(

 Create polynomial coefficient dataset from the loaded files.

 Reads and processes polynomial coefficient files to create a :class:`~tudatpy.dynamics.environment_setup.atmosphere.ComaPolyDataset`.
 This method is only available when the processor was constructed from polynomial coefficient files
 using the file path variant of :func:`~tudatpy.dynamics.environment_setup.atmosphere.coma_model_file_processor`.

 Returns
 -------
 ComaPolyDataset
     Dataset containing polynomial coefficients for all loaded files.

 Raises
 ------
 RuntimeError
     If the processor was constructed from Stokes coefficient CSV files instead of polynomial files.


 Examples
 --------
 .. code-block:: python

   # Create processor from polynomial files
   processor = environment_setup.atmosphere.coma_model_file_processor(
       ["h2o_poly_epoch1.txt", "h2o_poly_epoch2.txt"])

   # Create polynomial dataset
   poly_dataset = processor.create_poly_coef_dataset()

   # Use in coma atmosphere model
   coma_settings = environment_setup.atmosphere.coma_model(
       poly_data=poly_dataset,
       molecular_weight=0.018015)


      )doc")
        .def("create_sh_dataset",
             py::overload_cast<>(&tss::ComaModelFileProcessor::createSHDataset, py::const_),
             R"doc(

 Create Stokes coefficient dataset from preloaded CSV files (parameterless version).

 This method is only available when the processor was constructed from Stokes coefficient CSV files
 using :func:`~tudatpy.dynamics.environment_setup.atmosphere.coma_model_file_processor_from_sh_files`.
 It returns the preloaded Stokes dataset directly.

 Returns
 -------
 ComaStokesDataset
     Dataset containing preloaded Stokes coefficients from CSV files.

 Raises
 ------
 RuntimeError
     If the processor was constructed from polynomial coefficient files. Use the parameterized version
     :meth:`create_sh_dataset(radii_m, sol_longitudes_deg, ...)` for polynomial files instead.

 Examples
 --------
 .. code-block:: python

   # Create processor from Stokes CSV files
   processor = environment_setup.atmosphere.coma_model_file_processor_from_sh_files(
       input_dir="stokes_data",
       prefix="stokes")

   # Load Stokes dataset (no parameters needed)
   stokes_dataset = processor.create_sh_dataset()

   # Use in coma atmosphere model
   coma_settings = environment_setup.atmosphere.coma_model(
       stokes_data=stokes_dataset,
       molecular_weight=0.018015)

      )doc")
        .def("create_sh_dataset",
             py::overload_cast<const std::vector<double>&, const std::vector<double>&, const int, const int>(
                 &tss::ComaModelFileProcessor::createSHDataset, py::const_),
             py::arg("radii_m"),
             py::arg("sol_longitudes_deg"),
             py::arg("requested_max_degree") = -1,
             py::arg("requested_max_order") = -1,
             R"doc(

 Create Stokes coefficient dataset by transforming polynomial coefficients (parameterized version).

 This method is only available when the processor was constructed from polynomial coefficient files
 using :func:`~tudatpy.dynamics.environment_setup.atmosphere.coma_model_file_processor`. It transforms
 polynomial coefficients to Stokes coefficients by evaluating the spherical harmonic expansion at the
 specified grid of radii and solar longitudes.

 Parameters
 ----------
 radii_m : list[float]
     Vector of radii at which to evaluate Stokes coefficients [m].

 sol_longitudes_deg : list[float]
     Vector of solar longitudes at which to evaluate Stokes coefficients [degrees].

 requested_max_degree : int, default = -1
     Maximum spherical harmonic degree to include. Set to -1 to use maximum available.

 requested_max_order : int, default = -1
     Maximum spherical harmonic order to include. Set to -1 to use maximum available.

 Returns
 -------
 ComaStokesDataset
     Dataset containing Stokes coefficients.

 Raises
 ------
 RuntimeError
     If the processor was constructed from Stokes CSV files. Use the parameterless version
     :meth:`create_sh_dataset()` for Stokes CSV files instead.

 Examples
 --------
 .. code-block:: python

   # Create processor from polynomial files
   processor = environment_setup.atmosphere.coma_model_file_processor(
       ["h2o_poly.txt"])

   # Transform to Stokes at specific radii and solar longitudes
   stokes_dataset = processor.create_sh_dataset(
       radii_m=[1000.0, 2000.0, 5000.0, 10000.0, 20000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0],
       requested_max_degree=10,
       requested_max_order=10)

   # Use in coma atmosphere model
   coma_settings = environment_setup.atmosphere.coma_model(
       stokes_data=stokes_dataset,
       molecular_weight=0.018015)

      )doc")
        .def("create_sh_files",
             &tss::ComaModelFileProcessor::createSHFiles,
             py::arg("output_dir"),
             py::arg("radii_m"),
             py::arg("sol_longitudes_deg"),
             py::arg("requested_max_degree") = -1,
             py::arg("requested_max_order") = -1,
             R"doc(

 Create and save Stokes coefficient CSV files from polynomial coefficients.

 Transforms polynomial coefficients to Stokes coefficients and saves them as CSV files in the specified
 output directory. This is useful for pre-computing Stokes coefficients to avoid repeated transformations
 during multiple simulation runs. The generated CSV files can later be loaded using a processor created
 with the directory-based variant of :func:`~tudatpy.dynamics.environment_setup.atmosphere.coma_model_file_processor`.

 This method is only available when the processor was constructed from polynomial coefficient files.


 Parameters
 ----------
 output_dir : str
     Directory path where Stokes coefficient CSV files will be saved. The directory will be created
     if it does not exist.

 radii_m : list[float]
     Vector of radii at which to evaluate Stokes coefficients [m].

 sol_longitudes_deg : list[float]
     Vector of solar longitudes at which to evaluate Stokes coefficients [degrees].

 requested_max_degree : int, default = -1
     Maximum spherical harmonic degree to include. Set to -1 to automatically use the maximum
     degree available in the polynomial data.

 requested_max_order : int, default = -1
     Maximum spherical harmonic order to include. Set to -1 to automatically use the maximum
     order available in the polynomial data.

 Raises
 ------
 RuntimeError
     If the processor was constructed from Stokes CSV files instead of polynomial files.


 Examples
 --------
 .. code-block:: python

   # Create processor from polynomial files
   processor = environment_setup.atmosphere.coma_model_file_processor(
       ["h2o_poly_epoch1.txt", "h2o_poly_epoch2.txt"])

   # Pre-compute and save Stokes coefficients to CSV files
   processor.create_sh_files(
       output_dir="stokes_output",
       radii_m=[1000.0, 2000.0, 5000.0, 10000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0],
       requested_max_degree=10,
       requested_max_order=10)

   # Later, load from the saved CSV files
   processor_from_csv = environment_setup.atmosphere.coma_model_file_processor(
       input_dir="stokes_output",
       prefix="stokes")
   stokes_dataset = processor_from_csv.create_sh_dataset(
       radii_m=[],  # Ignored when loading from files
       sol_longitudes_deg=[])


      )doc");

    m.def(
            "coma_model_file_processor",
            py::overload_cast<const std::vector<std::string>&>(&tss::comaModelFileProcessorFromPolyFiles),
            py::arg( "file_paths" ),
            R"doc(

 Function for creating coma model file processor from polynomial coefficient files.

 Creates a :class:`~tudatpy.dynamics.environment_setup.atmosphere.ComaModelFileProcessor` that loads
 and processes polynomial coefficient files. These files contain spherical harmonic coefficients in
 polynomial form that describe gas density distributions in a cometary coma.

 The processor can create polynomial datasets directly, or transform them to Stokes coefficients at
 specified radii and solar longitudes. Multiple files can be provided to cover different time periods,
 enabling time-dependent coma modeling.


 Parameters
 ----------
 file_paths : list[str]
     List of file paths to polynomial coefficient files. Each file should contain polynomial coefficients
     for spherical harmonic expansion of gas density. Files may cover different time periods.

 Returns
 -------
 ComaModelFileProcessor
     Processor configured for polynomial coefficient files, capable of creating both polynomial and
     Stokes datasets.

 Raises
 ------
 ValueError
     If file_paths is an empty list.
 RuntimeError
     If any of the specified files do not exist or cannot be opened.


 Examples
 --------
 Create processor and use polynomial dataset directly:

 .. code-block:: python

   # Define paths to polynomial coefficient files
   poly_files = [
       "coma_data/h2o_poly_epoch1.txt",
       "coma_data/h2o_poly_epoch2.txt"
   ]

   # Create file processor
   processor = environment_setup.atmosphere.coma_model_file_processor(poly_files)

   # Create polynomial dataset
   poly_dataset = processor.create_poly_coef_dataset()

   # Create coma atmosphere settings
   coma_settings = environment_setup.atmosphere.coma_model(
       poly_data=poly_dataset,
       molecular_weight=0.018015,  # H2O in kg/mol
       max_degree=10,
       max_order=10)

   # Apply to body
   body_settings.get("67P").atmosphere_settings = coma_settings

 Transform to Stokes coefficients:

 .. code-block:: python

   # Create processor
   processor = environment_setup.atmosphere.coma_model_file_processor(poly_files)

   # Transform to Stokes coefficients at specific radii and solar longitudes
   stokes_dataset = processor.create_sh_dataset(
       radii_m=[1000.0, 2000.0, 5000.0, 10000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0],
       requested_max_degree=10,
       requested_max_order=10)

   # Create coma atmosphere settings
   coma_settings = environment_setup.atmosphere.coma_model(
       stokes_data=stokes_dataset,
       molecular_weight=0.018015)


      )doc"
            );

    m.def(
            "coma_model_file_processor",
            py::overload_cast<const std::string&, const std::string&>(&tss::comaModelFileProcessorFromSHFiles),
            py::arg( "input_dir" ),
            py::arg( "prefix" ) = "stokes",
            R"doc(

 Function for creating coma model file processor from Stokes coefficient CSV files.

 Creates a :class:`~tudatpy.dynamics.environment_setup.atmosphere.ComaModelFileProcessor` that loads
 pre-computed Stokes (spherical harmonic) coefficients from CSV files. These files typically contain
 Stokes coefficients that were previously generated and saved using the ``create_sh_files()`` method
 of a processor created from polynomial files.

 This variant is useful when you want to avoid repeated transformation of polynomial coefficients to
 Stokes coefficients across multiple simulation runs. The CSV files contain pre-evaluated coefficients
 at a fixed grid of radii and solar longitudes.


 Parameters
 ----------
 input_dir : str
     Directory path containing the Stokes coefficient CSV files to load. All CSV files in this
     directory with the specified prefix will be loaded.

 prefix : str, default = "stokes"
     File name prefix for the CSV files. Files are expected to be named as ``{prefix}_0.csv``,
     ``{prefix}_1.csv``, etc.

 Returns
 -------
 ComaModelFileProcessor
     Processor configured for Stokes coefficient CSV files, capable of loading pre-computed
     Stokes datasets.

 Raises
 ------
 RuntimeError
     If the directory does not exist, is not a directory, or contains no matching CSV files.


 Examples
 --------
 Load from previously saved Stokes CSV files:

 .. code-block:: python

   # Create processor from existing Stokes CSV files
   processor = environment_setup.atmosphere.coma_model_file_processor(
       input_dir="coma_data/stokes_precomputed",
       prefix="stokes")

   # Load Stokes dataset (radii and longitudes are read from files)
   stokes_dataset = processor.create_sh_dataset(
       radii_m=[],  # Ignored when loading from CSV files
       sol_longitudes_deg=[])  # Ignored when loading from CSV files

   # Create coma atmosphere settings
   coma_settings = environment_setup.atmosphere.coma_model(
       stokes_data=stokes_dataset,
       molecular_weight=0.018015)  # H2O in kg/mol

   # Apply to body
   body_settings.get("67P").atmosphere_settings = coma_settings

 Complete workflow from polynomial to saved Stokes to loading:

 .. code-block:: python

   # Step 1: Create and save Stokes coefficients from polynomial files
   poly_processor = environment_setup.atmosphere.coma_model_file_processor(
       ["h2o_poly.txt"])
   poly_processor.create_sh_files(
       output_dir="stokes_saved",
       radii_m=[1000.0, 2000.0, 5000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0])

   # Step 2: Later, load from saved Stokes files
   stokes_processor = environment_setup.atmosphere.coma_model_file_processor(
       input_dir="stokes_saved",
       prefix="stokes")
   stokes_dataset = stokes_processor.create_sh_dataset(
       radii_m=[],
       sol_longitudes_deg=[])

   # Step 3: Use in simulation
   coma_settings = environment_setup.atmosphere.coma_model(
       stokes_data=stokes_dataset,
       molecular_weight=0.018015)


      )doc"
            );

    // === Coma wind processing: file processor ===
    py::class_< tss::ComaWindModelFileProcessor, std::shared_ptr< tss::ComaWindModelFileProcessor > >( m,
                                               "ComaWindModelFileProcessor",
                                               R"doc(
        Processor for creating wind model datasets from three component file sources.

        This class manages the creation of ComaWindDatasetCollection from three sets of files
        (one for each spatial component: x, y, z). It provides a simplified interface for
        wind model setup by handling all three components together.
        )doc" )
        .def("create_poly_coef_datasets",
             &tss::ComaWindModelFileProcessor::createPolyCoefDatasets,
             R"doc(

 Create polynomial coefficient dataset collection for all three wind components.

 Reads and processes polynomial coefficient files for x, y, and z wind velocity components to create
 a :class:`~tudatpy.dynamics.environment_setup.atmosphere.ComaWindDatasetCollection`. This method is
 only available when the processor was constructed from polynomial coefficient files.

 Returns
 -------
 ComaWindDatasetCollection
     Collection containing x, y, z polynomial datasets for wind velocity components.

 Raises
 ------
 RuntimeError
     If processor was constructed from Stokes coefficient CSV files instead of polynomial files.


 Examples
 --------
 .. code-block:: python

   # Define paths to polynomial coefficient files for each component
   x_files = ["wind_x_epoch1.txt", "wind_x_epoch2.txt"]
   y_files = ["wind_y_epoch1.txt", "wind_y_epoch2.txt"]
   z_files = ["wind_z_epoch1.txt", "wind_z_epoch2.txt"]

   # Create wind file processor
   wind_processor = environment_setup.atmosphere.coma_wind_file_processor(
       x_files, y_files, z_files)

   # Create polynomial dataset collection
   poly_datasets = wind_processor.create_poly_coef_datasets()


      )doc")
        .def("create_sh_datasets",
             py::overload_cast<>(&tss::ComaWindModelFileProcessor::createSHDatasets, py::const_),
             R"doc(

 Create Stokes coefficient dataset collection from preloaded CSV files (parameterless version).

 This method is only available when the processor was constructed from Stokes coefficient CSV files
 using :func:`~tudatpy.dynamics.environment_setup.atmosphere.coma_wind_file_processor_from_sh_files`.
 It returns the preloaded Stokes datasets for all three wind components (x, y, z) directly.

 Returns
 -------
 ComaWindDatasetCollection
     Collection containing preloaded x, y, z Stokes datasets for wind velocity components.

 Raises
 ------
 RuntimeError
     If the processor was constructed from polynomial coefficient files. Use the parameterized version
     :meth:`create_sh_datasets(radii_m, sol_longitudes_deg, ...)` for polynomial files instead.

 Examples
 --------
 .. code-block:: python

   # Create wind processor from Stokes CSV files
   wind_processor = environment_setup.atmosphere.coma_wind_file_processor_from_sh_files(
       x_input_dir="wind_x_stokes",
       y_input_dir="wind_y_stokes",
       z_input_dir="wind_z_stokes",
       prefix="stokes")

   # Load Stokes datasets (no parameters needed)
   stokes_datasets = wind_processor.create_sh_datasets()

   # Use in coma wind model
   coma_wind = environment_setup.atmosphere.coma_wind_model(
       dataset_collection=stokes_datasets,
       requested_max_degree=10,
       requested_max_order=10)

      )doc")
        .def("create_sh_datasets",
             py::overload_cast<const std::vector<double>&, const std::vector<double>&, const int, const int>(
                 &tss::ComaWindModelFileProcessor::createSHDatasets, py::const_),
             py::arg("radii_m"),
             py::arg("sol_longitudes_deg"),
             py::arg("requested_max_degree") = -1,
             py::arg("requested_max_order") = -1,
             R"doc(

 Create Stokes coefficient dataset collection by transforming polynomial coefficients (parameterized version).

 This method is only available when the processor was constructed from polynomial coefficient files
 using :func:`~tudatpy.dynamics.environment_setup.atmosphere.coma_wind_file_processor`. It transforms
 polynomial coefficients to Stokes coefficients for all three components (x, y, z) by evaluating
 at the specified grid of radii and solar longitudes.

 Parameters
 ----------
 radii_m : list[float]
     Vector of radii at which to evaluate Stokes coefficients [m].

 sol_longitudes_deg : list[float]
     Vector of solar longitudes at which to evaluate Stokes coefficients [degrees].

 requested_max_degree : int, default = -1
     Maximum spherical harmonic degree to include (-1 for automatic determination).

 requested_max_order : int, default = -1
     Maximum spherical harmonic order to include (-1 for automatic determination).

 Returns
 -------
 ComaWindDatasetCollection
     Collection containing x, y, z Stokes datasets for wind velocity components.

 Raises
 ------
 RuntimeError
     If the processor was constructed from Stokes CSV files. Use the parameterless version
     :meth:`create_sh_datasets()` for Stokes CSV files instead.

 Examples
 --------
 .. code-block:: python

   # Create wind processor from polynomial files
   wind_processor = environment_setup.atmosphere.coma_wind_file_processor(
       x_file_paths=["wind_x.txt"],
       y_file_paths=["wind_y.txt"],
       z_file_paths=["wind_z.txt"])

   # Transform to Stokes coefficients at specific radii and solar longitudes
   stokes_datasets = wind_processor.create_sh_datasets(
       radii_m=[1000.0, 2000.0, 5000.0, 10000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0],
       requested_max_degree=10,
       requested_max_order=10)

   # Use in coma wind model
   coma_wind = environment_setup.atmosphere.coma_wind_model(
       dataset_collection=stokes_datasets,
       requested_max_degree=10,
       requested_max_order=10)

      )doc")
        .def("create_sh_files",
             &tss::ComaWindModelFileProcessor::createSHFiles,
             py::arg("x_output_dir"),
             py::arg("y_output_dir"),
             py::arg("z_output_dir"),
             py::arg("radii_m"),
             py::arg("sol_longitudes_deg"),
             py::arg("requested_max_degree") = -1,
             py::arg("requested_max_order") = -1,
             R"doc(

 Create and save Stokes coefficient CSV files for all three wind components.

 Transforms polynomial coefficients to Stokes coefficients for x, y, and z wind velocity components
 and saves them as CSV files in separate output directories. This is useful for pre-computing Stokes
 coefficients to avoid repeated transformations across multiple simulation runs.

 This method is only available when the processor was constructed from polynomial coefficient files.


 Parameters
 ----------
 x_output_dir : str
     Directory path where x-component Stokes coefficient CSV files will be saved. The directory
     will be created if it does not exist.

 y_output_dir : str
     Directory path where y-component Stokes coefficient CSV files will be saved. The directory
     will be created if it does not exist.

 z_output_dir : str
     Directory path where z-component Stokes coefficient CSV files will be saved. The directory
     will be created if it does not exist.

 radii_m : list[float]
     Vector of radii at which to evaluate Stokes coefficients [m].

 sol_longitudes_deg : list[float]
     Vector of solar longitudes at which to evaluate Stokes coefficients [degrees].

 requested_max_degree : int, default = -1
     Maximum spherical harmonic degree to include (-1 for automatic determination from data).

 requested_max_order : int, default = -1
     Maximum spherical harmonic order to include (-1 for automatic determination from data).

 Raises
 ------
 RuntimeError
     If processor was constructed from Stokes CSV files instead of polynomial files.


 Examples
 --------
 .. code-block:: python

   # Create wind processor from polynomial files
   wind_processor = environment_setup.atmosphere.coma_wind_file_processor(
       x_file_paths=["wind_x_epoch1.txt"],
       y_file_paths=["wind_y_epoch1.txt"],
       z_file_paths=["wind_z_epoch1.txt"])

   # Pre-compute and save Stokes coefficients for all components
   wind_processor.create_sh_files(
       x_output_dir="stokes_wind/x_component",
       y_output_dir="stokes_wind/y_component",
       z_output_dir="stokes_wind/z_component",
       radii_m=[1000.0, 2000.0, 5000.0, 10000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0],
       requested_max_degree=10,
       requested_max_order=10)

   # Later, load from the saved CSV files
   wind_processor_from_csv = environment_setup.atmosphere.coma_wind_file_processor(
       x_input_dir="stokes_wind/x_component",
       y_input_dir="stokes_wind/y_component",
       z_input_dir="stokes_wind/z_component",
       prefix="stokes")


      )doc")
        .def("is_poly_type", &tss::ComaWindModelFileProcessor::isPolyType,
             R"doc(Check if this processor works with polynomial coefficient files.)doc")
        .def("is_stokes_type", &tss::ComaWindModelFileProcessor::isStokesType,
             R"doc(Check if this processor works with Stokes coefficient files.)doc");

    m.def(
            "coma_wind_file_processor",
            py::overload_cast<const std::vector<std::string>&,
                              const std::vector<std::string>&,
                              const std::vector<std::string>&>(&tss::comaWindModelFileProcessorFromPolyFiles),
            py::arg( "x_file_paths" ),
            py::arg( "y_file_paths" ),
            py::arg( "z_file_paths" ),
            R"doc(

 Function for creating coma wind model file processor from polynomial coefficient files.

 Creates a :class:`~tudatpy.dynamics.environment_setup.atmosphere.ComaWindModelFileProcessor` that loads
 and processes polynomial coefficient files for all three wind velocity components (x, y, z). These files
 contain spherical harmonic coefficients in polynomial form that describe wind velocity distributions
 in a cometary coma.

 The processor can create polynomial dataset collections directly, or transform them to Stokes
 coefficients at specified radii and solar longitudes. Multiple files can be provided for each component
 to cover different time periods.


 Parameters
 ----------
 x_file_paths : list[str]
     List of file paths for x-component polynomial coefficients. Files may cover different time periods.

 y_file_paths : list[str]
     List of file paths for y-component polynomial coefficients. Files may cover different time periods.

 z_file_paths : list[str]
     List of file paths for z-component polynomial coefficients. Files may cover different time periods.

 Returns
 -------
 ComaWindModelFileProcessor
     Processor configured for polynomial coefficient files, capable of creating both polynomial and
     Stokes dataset collections.

 Raises
 ------
 ValueError
     If any of the file path lists is empty.
 RuntimeError
     If any of the specified files do not exist or cannot be opened.


 Examples
 --------
 Create wind model from polynomial files:

 .. code-block:: python

   # Define paths to polynomial coefficient files for each wind component
   x_files = ["wind_data/vx_epoch1.txt", "wind_data/vx_epoch2.txt"]
   y_files = ["wind_data/vy_epoch1.txt", "wind_data/vy_epoch2.txt"]
   z_files = ["wind_data/vz_epoch1.txt", "wind_data/vz_epoch2.txt"]

   # Create wind file processor
   wind_processor = environment_setup.atmosphere.coma_wind_file_processor(
       x_file_paths=x_files,
       y_file_paths=y_files,
       z_file_paths=z_files)

   # Transform to Stokes coefficients
   wind_datasets = wind_processor.create_sh_datasets(
       radii_m=[1000.0, 2000.0, 5000.0, 10000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0],
       requested_max_degree=10,
       requested_max_order=10)

   # Create coma wind model settings
   coma_wind = environment_setup.atmosphere.coma_wind_model(
       dataset_collection=wind_datasets,
       requested_max_degree=10,
       requested_max_order=10,
       associated_reference_frame=environment.AerodynamicsReferenceFrames.vertical_frame,
       include_corotation=True)

   # Apply to atmosphere
   body_settings.get("67P").atmosphere_settings.wind_settings = coma_wind

 Pre-compute and save Stokes coefficients:

 .. code-block:: python

   # Create processor
   wind_processor = environment_setup.atmosphere.coma_wind_file_processor(
       x_file_paths=x_files,
       y_file_paths=y_files,
       z_file_paths=z_files)

   # Save Stokes coefficients to CSV files for later use
   wind_processor.create_sh_files(
       x_output_dir="stokes_wind/x",
       y_output_dir="stokes_wind/y",
       z_output_dir="stokes_wind/z",
       radii_m=[1000.0, 2000.0, 5000.0, 10000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0],
       requested_max_degree=10,
       requested_max_order=10)


      )doc"
            );

    m.def(
            "coma_wind_file_processor",
            py::overload_cast<const std::string&,
                              const std::string&,
                              const std::string&,
                              const std::string&>(&tss::comaWindModelFileProcessorFromSHFiles),
            py::arg( "x_input_dir" ),
            py::arg( "y_input_dir" ),
            py::arg( "z_input_dir" ),
            py::arg( "prefix" ) = "stokes",
            R"doc(

 Function for creating coma wind model file processor from Stokes coefficient CSV files.

 Creates a :class:`~tudatpy.dynamics.environment_setup.atmosphere.ComaWindModelFileProcessor` that loads
 pre-computed Stokes (spherical harmonic) coefficients from CSV files for all three wind velocity
 components (x, y, z). These files typically contain Stokes coefficients that were previously generated
 and saved using the ``create_sh_files()`` method of a processor created from polynomial files.

 This variant is useful when you want to avoid repeated transformation of polynomial coefficients to
 Stokes coefficients across multiple simulation runs. The CSV files contain pre-evaluated coefficients
 at a fixed grid of radii and solar longitudes.


 Parameters
 ----------
 x_input_dir : str
     Directory path containing x-component Stokes coefficient CSV files. All CSV files in this
     directory with the specified prefix will be loaded.

 y_input_dir : str
     Directory path containing y-component Stokes coefficient CSV files. All CSV files in this
     directory with the specified prefix will be loaded.

 z_input_dir : str
     Directory path containing z-component Stokes coefficient CSV files. All CSV files in this
     directory with the specified prefix will be loaded.

 prefix : str, default = "stokes"
     File name prefix for the CSV files. Files are expected to be named as ``{prefix}_0.csv``,
     ``{prefix}_1.csv``, etc.

 Returns
 -------
 ComaWindModelFileProcessor
     Processor configured for Stokes coefficient CSV files, capable of loading pre-computed
     Stokes dataset collections.

 Raises
 ------
 RuntimeError
     If any of the directories do not exist, are not directories, or contain no matching CSV files.


 Examples
 --------
 Load from previously saved Stokes CSV files:

 .. code-block:: python

   # Create processor from existing Stokes CSV files
   wind_processor = environment_setup.atmosphere.coma_wind_file_processor(
       x_input_dir="wind_stokes/x_component",
       y_input_dir="wind_stokes/y_component",
       z_input_dir="wind_stokes/z_component",
       prefix="stokes")

   # Load Stokes datasets (radii and longitudes are read from files)
   wind_datasets = wind_processor.create_sh_datasets(
       radii_m=[],  # Ignored when loading from CSV files
       sol_longitudes_deg=[])  # Ignored when loading from CSV files

   # Create coma wind model settings
   coma_wind = environment_setup.atmosphere.coma_wind_model(
       dataset_collection=wind_datasets,
       associated_reference_frame=environment.AerodynamicsReferenceFrames.vertical_frame)

   # Apply to atmosphere
   body_settings.get("67P").atmosphere_settings.wind_settings = coma_wind

 Complete workflow from polynomial to saved Stokes to loading:

 .. code-block:: python

   # Step 1: Create and save Stokes coefficients from polynomial files
   poly_processor = environment_setup.atmosphere.coma_wind_file_processor(
       x_file_paths=["vx.txt"],
       y_file_paths=["vy.txt"],
       z_file_paths=["vz.txt"])
   poly_processor.create_sh_files(
       x_output_dir="wind_stokes/x",
       y_output_dir="wind_stokes/y",
       z_output_dir="wind_stokes/z",
       radii_m=[1000.0, 2000.0, 5000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0])

   # Step 2: Later, load from saved Stokes files
   stokes_processor = environment_setup.atmosphere.coma_wind_file_processor(
       x_input_dir="wind_stokes/x",
       y_input_dir="wind_stokes/y",
       z_input_dir="wind_stokes/z",
       prefix="stokes")
   wind_datasets = stokes_processor.create_sh_datasets(
       radii_m=[],
       sol_longitudes_deg=[])

   # Step 3: Use in simulation
   coma_wind = environment_setup.atmosphere.coma_wind_model(
       dataset_collection=wind_datasets)
   body_settings.get("67P").atmosphere_settings.wind_settings = coma_wind


      )doc"
            );



    m.def("mars_dtm",
          &tss::marsDtmAtmosphereSettings,
          R"doc(No documentation found.)doc" );
}

}  // namespace atmosphere
}  // namespace environment_setup
}  // namespace dynamics
}  // namespace tudatpy
