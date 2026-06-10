/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif
#include "expose_atmosphere.h"

#include <tudat/astro/aerodynamics/nrlmsise00Atmosphere.h>
#include <tudat/astro/aerodynamics/nrlmsise00InputFunctions.h>
#include <tudat/astro/reference_frames/referenceFrameTransformations.h>
#include <tudat/astro/aerodynamics/mcdAtmosphereModel.h>
#include <tudat/simulation/environment_setup/createAtmosphereModel.h>
#include <tudat/paths.hpp>

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
    std::string atmosphereTableFile = paths::getAtmosphereTablesPath( ) + "/USSA1976Until100kmPer100mUntil1000kmPer1000m.dat";
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

void expose_atmosphere_setup( py::module& m )
{
    // NRLMSISE00
    py::class_< ta::NRLMSISE00Input, std::shared_ptr< ta::NRLMSISE00Input > >( m,
                                                                               "NRLMSISE00Input",
                                                                               R"doc(
Input for computation of NRLMSISE00 atmospheric conditions at current time and position.

The values in this class may be recomputed every time step to reflect changing atmospheric conditions.
)doc" )
            .def( py::init< int, int, double, double, double, double, double, std::vector< double >, std::vector< int > >( ),
                  py::arg( "year" ) = 0,
                  py::arg( "day_of_year" ) = 0,
                  py::arg( "seconds_of_day" ) = 0.0,
                  py::arg( "local_solar_time" ) = 0.0,
                  py::arg( "f107" ) = 0.0,
                  py::arg( "f107a" ) = 0.0,
                  py::arg( "ap_daily" ) = 0.0,
                  py::arg( "ap_vector" ) = std::vector< double >( 7, 0.0 ),
                  py::arg( "switches" ) = std::vector< int >( ) );

    py::class_< ta::NRLMSISE00Atmosphere, std::shared_ptr< ta::NRLMSISE00Atmosphere > >( m,
                                                                                         "NRLMSISE00Atmosphere",
                                                                                         R"doc(NRLMSISE00 atmosphere model.

                         This class uses the NRLMSISE00 model to compute the atmospheric density and temperature. The GTD7 function is used: Neutral Atmosphere Empirical Model from the surface to the lower exosphere.

                         Currently, the ideal gas law is used to compute the speed of sound and the specific heat ratio is assumed to be constant and equal to 1.4.

                         :param solar_activity_data: Solar activity data for a range of epochs as produced by tudatpy.data.read_solar_activity_data.
                         )doc" )
            .def( py::init< const std::map< double, std::shared_ptr< tio::solar_activity::SolarActivityData > >,
                            const bool,
                            const bool,
                            const bool >( ),
                  py::arg( "solar_activity_data" ),
                  py::arg( "use_ideal_gas_law" ) = true,
                  py::arg( "use_storm_conditions" ) = false,
                  py::arg( "use_anomalous_oxygen" ) = true )
            .def( "set_use_geodetic_latitude", &ta::NRLMSISE00Atmosphere::setUseGeodeticLatitude )
            .def( "get_use_geodetic_latitude", &ta::NRLMSISE00Atmosphere::getUseGeodeticLatitude )
            .def( "set_use_utc", &ta::NRLMSISE00Atmosphere::setUseUtc )
            .def( "get_use_utc", &ta::NRLMSISE00Atmosphere::getUseUtc )
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
            .value( "tabulated_density", tss::AtmosphereDependentVariables::density_dependent_atmosphere )
            .value( "tabulated_pressure", tss::AtmosphereDependentVariables::pressure_dependent_atmosphere )
            .value( "tabulated_temperature", tss::AtmosphereDependentVariables::temperature_dependent_atmosphere )
            .value( "tabulated_gas_constant", tss::AtmosphereDependentVariables::gas_constant_dependent_atmosphere )
            .value( "tabulated_specific_heat_ratio", tss::AtmosphereDependentVariables::specific_heat_ratio_dependent_atmosphere )
            .value( "tabulated_molar_mass", tss::AtmosphereDependentVariables::molar_mass_dependent_atmosphere )
            .export_values( );

    /////////////////////////////////////////////////////////////////////////////
    py::class_< tss::WindModelSettings, std::shared_ptr< tss::WindModelSettings > >( m,
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

    py::class_< tss::CustomWindModelSettings, std::shared_ptr< tss::CustomWindModelSettings >, tss::WindModelSettings >(
            m, "CustomWindModelSettings", R"doc(No documentation found.)doc" );

    py::class_< tss::AtmosphereSettings, std::shared_ptr< tss::AtmosphereSettings > >( m,
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

    py::class_< tss::ExponentialAtmosphereSettings, std::shared_ptr< tss::ExponentialAtmosphereSettings >, tss::AtmosphereSettings >(
            m,
            "ExponentialAtmosphereSettings",
            R"doc(

         Class for providing settings for exponential atmosphere model.

         `AtmosphereSettings` derived class for a defining the settings of an exponential atmosphere model.




      )doc" );

    py::class_< tss::CustomConstantTemperatureAtmosphereSettings,
                std::shared_ptr< tss::CustomConstantTemperatureAtmosphereSettings >,
                tss::AtmosphereSettings >( m, "CustomConstantTemperatureAtmosphereSettings", R"doc(No documentation found.)doc" );

    py::class_< tss::CustomNumberDensityAtmosphereSettings,
                std::shared_ptr< tss::CustomNumberDensityAtmosphereSettings >,
                tss::AtmosphereSettings >( m, "CustomNumberDensityAtmosphereSettings", R"doc(No documentation found.)doc" );

    py::class_< tss::ScaledAtmosphereSettings, std::shared_ptr< tss::ScaledAtmosphereSettings >, tss::AtmosphereSettings >(
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
 wind velocity components in a modified vertical frame. The wind model uses spherical harmonic
 expansion to compute wind velocities as a function of position.

 .. important::
     **Data fitting requirement**: The polynomial/Stokes coefficients for the wind model must be
     fitted from **raw (untransformed)** wind velocity values in m/s. Unlike the coma density model
     (which uses log2-transformed data), the wind model operates directly on the actual velocity values
     without any logarithmic transformation.

 The wind velocity components are defined in a **modified vertical frame**:

 * **X-axis**: Meridional direction (in the meridian plane, pointing towards the North, aligned with central-body-fixed Z-axis direction)
 * **Y-axis**: Zonal direction (completes the right-handed frame, pointing towards the West)
 * **Z-axis**: Radial direction pointing **OUTWARD** from the comet nucleus center (away from origin)

 .. warning::
     The Z-axis direction is **OPPOSITE** to the standard Tudat vertical frame convention, where
     Z points inward along the gravity vector. For the coma wind model, positive Z points radially
     outward. This is critical for correct sign conventions when preparing input data.


 Parameters
 ----------
 dataset_collection : ComaWindDatasetCollection
     Collection containing wind component datasets in the modified vertical frame (either polynomial or Stokes coefficients).

 requested_max_degree : int, default = -1
     Maximum spherical harmonic degree to use (-1 for automatic determination from data).

 requested_max_order : int, default = -1
     Maximum spherical harmonic order to use (-1 for automatic determination from data).

 associated_reference_frame : dynamics.environment.AerodynamicsReferenceFrames, default = AerodynamicsReferenceFrames.vertical_frame
     Reference frame in which the wind velocity is defined. For coma wind model, this uses a modified
     vertical frame with Z-axis pointing radially outward (away from nucleus), opposite to standard vertical frame.

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
   wind_processor = data.coma_model.coma_wind_file_processor(
       x_file_paths, y_file_paths, z_file_paths)

   # Create dataset collection with Stokes coefficients
   wind_datasets = wind_processor.create_coma_stokes_dataset(
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
           py::overload_cast< const std::string& >( &tss::exponentialAtmosphereSettings ),
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
           py::overload_cast< const double, const double, const double, const double, const double >( &tss::exponentialAtmosphereSettings ),
           py::arg( "scale_height" ),
           py::arg( "surface_density" ),
           py::arg( "constant_temperature" ) = 288.15,
           py::arg( "specific_gas_constant" ) = tudat::physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
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
           py::arg( "space_weather_file" ) = tudat::paths::getSpaceWeatherDataPath( ) + "/sw19571001.txt",
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
                   { tss::density_dependent_atmosphere, tss::pressure_dependent_atmosphere, tss::temperature_dependent_atmosphere } ),
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
           py::overload_cast< const std::function< double( const double ) >, const double, const double, const double >(
                   &tss::customConstantTemperatureAtmosphereSettings ),
           py::arg( "density_function" ),
           py::arg( "constant_temperature" ),
           py::arg( "specific_gas_constant" ) = tudat::physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
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
           py::overload_cast< const std::function< double( const double, const double, const double, const double ) >,
                              const double,
                              const double,
                              const double >( &tss::customConstantTemperatureAtmosphereSettings ),
           py::arg( "density_function" ),
           py::arg( "constant_temperature" ),
           py::arg( "specific_gas_constant" ) = tudat::physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
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

    m.def( "custom_number_density",
           py::overload_cast< const std::function< double( const double, const double, const double, const double ) >,
                              const double,
                              const double,
                              const double >( &tss::customNumberDensityAtmosphereSettings ),
           py::arg( "number_density_function" ),
           py::arg( "molar_mass" ),
           py::arg( "constant_temperature" ) = TUDAT_NAN,
           py::arg( "ratio_of_specific_heats" ) = 1.4,
           R"doc(

 Function for creating atmospheric model settings from a custom number-density profile.

 The supplied function returns total number density in m^-3 as a function of altitude,
 longitude, latitude and time. Mass density is computed internally from the molar mass
 using Avogadro's constant.

 Parameters
 ----------
 number_density_function : callable[[float, float, float, float], float]
     Function to retrieve the total number density at the current altitude, longitude,
     latitude and time.
 molar_mass : float
     Molar mass of the atmospheric species in kg/mol.
 constant_temperature : float, optional
     Constant temperature used only for pressure and speed-of-sound queries.
 ratio_of_specific_heats : float, default = 1.4
     Ratio of specific heats used only for speed-of-sound queries.
 Returns
 -------
 CustomNumberDensityAtmosphereSettings
     Settings for a custom number-density-driven atmosphere.

     )doc" );

    m.def( "scaled_by_function",
           py::overload_cast< const std::shared_ptr< tss::AtmosphereSettings >, const std::function< double( const double ) >, const bool >(
                   &tss::scaledAtmosphereSettings ),
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
 density_scaling_function : callable[[float], float]
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
           py::overload_cast< const std::shared_ptr< tss::AtmosphereSettings >, const double, const bool >(
                   &tss::scaledAtmosphereSettings ),
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
            "coma_model_from_poly_data",
            py::overload_cast<
                const tss::ComaPolyDataset&, double, int, int, bool >( &tss::comaSettings ),
            py::arg( "poly_data" ),
            py::arg( "molecular_weight" ),
            py::arg( "max_degree" ) = -1,
            py::arg( "max_order" ) = -1,
            py::arg( "is_log2" ) = true,
            R"doc(

 Function for creating coma atmosphere model settings from polynomial coefficients.

 Function for settings object, defining a coma atmosphere model based on spherical harmonic expansion
 of gas density data. The coma model is designed for modeling cometary atmospheres (comae) where
 gas density varies with position and time. This variant uses polynomial coefficient data that
 describes the spatial distribution of gas density.

 The density is computed using spherical harmonic expansion, allowing efficient representation of
 complex 3D density distributions around the nucleus. The model supports time-dependent density
 variations through multiple data files covering different time periods.

 .. note::
     **Data fitting**: By default (``is_log2=True``), the model assumes that the polynomial coefficients
     were fitted from **log2-transformed** number density data (i.e., log2(n) where n is the number
     density in m^-3) and internally applies the inverse transformation (2^x). If your coefficients
     were fitted to raw (non-log2) number density data, set ``is_log2=False`` to skip the
     back-transformation.


 Parameters
 ----------
 poly_data : ComaPolyDataset
     Polynomial coefficient dataset containing spherical harmonic coefficients for gas density
     distribution. Create using :func:`~tudatpy.data.coma_model.coma_model_file_processor`.

 molecular_weight : float
     Molecular weight (molar mass) of the gas species [kg/mol]. For water vapor (H2O), use 0.018015 kg/mol.

 max_degree : int, default = -1
     Maximum spherical harmonic degree to use in density calculations. Set to -1 to automatically
     use the maximum degree available in the dataset.

 max_order : int, default = -1
     Maximum spherical harmonic order to use in density calculations. Set to -1 to automatically
     use the maximum order available in the dataset.

 is_log2 : bool, default = True
     Whether the coefficients were fitted to log2-transformed number density data. If True (default),
     the model applies exp2 to convert back to actual number density. If False, the spherical
     harmonics output is used directly as number density.

 Returns
 -------
 ComaSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.ComaSettings` class
     configured for coma model. The returned object can be used to add a temperature model via
     the :meth:`~tudatpy.dynamics.environment_setup.atmosphere.ComaSettings.add_temperature_model` method.


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
   processor = data.coma_model.coma_model_file_processor(poly_file_paths)

   # Create polynomial dataset
   poly_dataset = processor.create_poly_coefficient_dataset()

   # Create coma atmosphere settings
   coma_settings = environment_setup.atmosphere.coma_model_from_poly_data(
       poly_data=poly_dataset,
       molecular_weight=0.018015,  # H2O molecular weight in kg/mol
       max_degree=10,
       max_order=10)

   # Optionally add temperature model
   # coma_settings.add_temperature_model(poly_data=temperature_poly_dataset)

   # Apply to body settings
   body_settings.get("67P").atmosphere_settings = coma_settings


    )doc"
                );

    m.def(
            "coma_model_from_stokes_data",
            py::overload_cast<
                const tss::ComaStokesDataset&, double, int, int, bool >( &tss::comaSettings ),
            py::arg( "stokes_data" ),
            py::arg( "molecular_weight" ),
            py::arg( "max_degree" ) = -1,
            py::arg( "max_order" ) = -1,
            py::arg( "is_log2" ) = true,
            R"doc(

 Function for creating coma atmosphere model settings from Stokes coefficients.

 Function for settings object, defining a coma atmosphere model based on spherical harmonic expansion
 of gas density data using precomputed Stokes coefficients. The coma model is designed for modeling
 cometary atmospheres (comae) where gas density varies with position and time. This variant uses
 precomputed Stokes coefficients (spherical harmonics) evaluated at specific radii and solar longitudes.

 Stokes coefficients provide a more direct representation of the spherical harmonic expansion compared
 to polynomial coefficients, offering faster evaluation during simulation. The coefficients
 are pre-evaluated at a grid of radii and solar longitudes, with interpolation used for intermediate values.

 .. note::
     **Data fitting**: By default (``is_log2=True``), the model assumes that the Stokes coefficients
     (or the polynomial coefficients they are derived from) were fitted from **log2-transformed**
     number density data (i.e., log2(n) where n is the number density in m^-3) and internally
     applies the inverse transformation (2^x). If your coefficients were fitted to raw (non-log2)
     number density data, set ``is_log2=False`` to skip the back-transformation.


 Parameters
 ----------
 stokes_data : ComaStokesDataset
     Precomputed Stokes coefficient dataset containing spherical harmonic coefficients evaluated at
     specific radii and solar longitudes. Create using :func:`~tudatpy.data.coma_model.coma_model_file_processor`
     or load from pre-existing Stokes coefficient CSV files.

 molecular_weight : float
     Molecular weight (molar mass) of the gas species [kg/mol]. For water vapor (H2O), use 0.018015 kg/mol.

 max_degree : int, default = -1
     Maximum spherical harmonic degree to use in density calculations. Set to -1 to automatically
     use the maximum degree available in the dataset.

 max_order : int, default = -1
     Maximum spherical harmonic order to use in density calculations. Set to -1 to automatically
     use the maximum order available in the dataset.

 is_log2 : bool, default = True
     Whether the coefficients were fitted to log2-transformed number density data. If True (default),
     the model applies exp2 to convert back to actual number density. If False, the spherical
     harmonics output is used directly as number density.

 Returns
 -------
 ComaSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.ComaSettings` class
     configured for coma model. The returned object can be used to add a temperature model via
     the :meth:`~tudatpy.dynamics.environment_setup.atmosphere.ComaSettings.add_temperature_model` method.


 Examples
 --------
 In this example, we create a coma atmosphere model by converting polynomial coefficients to Stokes coefficients:

 .. code-block:: python

   # Create file processor from polynomial coefficient files
   poly_file_paths = ["coma_data/poly_coeffs_epoch1.txt"]
   processor = data.coma_model.coma_model_file_processor(poly_file_paths)

   # Create Stokes dataset by evaluating at specific radii and solar longitudes
   stokes_dataset = processor.create_coma_stokes_dataset(
       radii_m=[1000.0, 2000.0, 5000.0, 10000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0],
       requested_max_degree=10,
       requested_max_order=10)

   # Create coma atmosphere settings
   coma_settings = environment_setup.atmosphere.coma_model_from_stokes_data(
       stokes_data=stokes_dataset,
       molecular_weight=0.018015,  # H2O molecular weight in kg/mol
       max_degree=10,
       max_order=10)

   # Optionally add temperature model
   # coma_settings.add_temperature_model(stokes_data=temperature_stokes_dataset)

   # Apply to body settings
   body_settings.get("67P").atmosphere_settings = coma_settings

 Alternatively, load from pre-existing Stokes coefficient CSV files:

 .. code-block:: python

   # Create file processor from existing Stokes CSV files
   processor = data.coma_model.coma_model_file_processor(
       input_dir="coma_data/stokes_files",
       prefix="stokes")

   # Create Stokes dataset (radii and longitudes are read from files)
   stokes_dataset = processor.create_coma_stokes_dataset(
       radii_m=[],  # Ignored when loading from files
       sol_longitudes_deg=[])

   # Create coma atmosphere settings
   coma_settings = environment_setup.atmosphere.coma_model_from_stokes_data(
       stokes_data=stokes_dataset,
       molecular_weight=0.018015)

   # Optionally add temperature model
   # coma_settings.add_temperature_model(stokes_data=temperature_stokes_dataset)

   # Apply to body settings
   body_settings.get("67P").atmosphere_settings = coma_settings


    )doc"
            );

    // === ComaSettings class exposure ===
    py::class_< tss::ComaSettings,
                std::shared_ptr< tss::ComaSettings >,
                tss::AtmosphereSettings >(
            m,
            "ComaSettings",
            R"doc(
Settings class for coma atmosphere models.

This class extends :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings`
to provide configuration for cometary coma atmosphere models. It supports both polynomial
and Stokes coefficient datasets for density modeling, and allows optional temperature
modeling via the :meth:`add_temperature_model` method.

.. note:: This class is typically created using the factory functions
          :func:`~tudatpy.dynamics.environment_setup.atmosphere.coma_model` rather than
          being instantiated directly.

Examples
--------
Create coma settings and add temperature model:

.. code-block:: python

  # Create coma atmosphere settings from polynomial data
  coma_settings = environment_setup.atmosphere.coma_model_from_poly_data(
      poly_data=density_poly_dataset,
      molecular_weight=0.018015)  # H2O in kg/mol

  # Add temperature model using polynomial data
  coma_settings.add_temperature_model(
      poly_data=temperature_poly_dataset,
      max_degree=10,
      max_order=10,
      gamma=1.33)

  # Apply to body
  body_settings.get("67P").atmosphere_settings = coma_settings

)doc" )
            .def( "add_temperature_model",
                  py::overload_cast< const tss::ComaPolyDataset&, const int, const int, const double >(
                          &tss::ComaSettings::addTemperatureModel ),
                  py::arg( "poly_data" ),
                  py::arg( "max_degree" ) = -1,
                  py::arg( "max_order" ) = -1,
                  py::arg( "gamma" ) = 1.33,
                  R"doc(
Add temperature model from polynomial coefficient data.

This method adds a temperature model to the coma atmosphere settings using polynomial
coefficient data. The temperature model uses spherical harmonic expansion to compute
temperature as a function of position and time.

.. note:: The temperature data type (polynomial or Stokes) must match the density data type.

Parameters
----------
poly_data : ComaPolyDataset
    Polynomial coefficient dataset for temperature distribution.

max_degree : int, default = -1
    Maximum spherical harmonic degree for temperature calculations. Set to -1 to use
    the maximum degree available in the dataset.

max_order : int, default = -1
    Maximum spherical harmonic order for temperature calculations. Set to -1 to use
    the maximum order available in the dataset.

gamma : float, default = 1.33
    Heat capacity ratio (gamma = Cp/Cv) for the gas species. Default value 1.33 is
    appropriate for water vapor.

Examples
--------
.. code-block:: python

  # Create coma settings with polynomial density data
  coma_settings = environment_setup.atmosphere.coma_model_from_poly_data(
      poly_data=density_poly_data,
      molecular_weight=0.018015)

  # Add temperature model with polynomial data
  coma_settings.add_temperature_model(
      poly_data=temperature_poly_data,
      max_degree=10,
      max_order=10,
      gamma=1.33)

)doc" )
            .def( "add_temperature_model",
                  py::overload_cast< const tss::ComaStokesDataset&, const int, const int, const double >(
                          &tss::ComaSettings::addTemperatureModel ),
                  py::arg( "stokes_data" ),
                  py::arg( "max_degree" ) = -1,
                  py::arg( "max_order" ) = -1,
                  py::arg( "gamma" ) = 1.33,
                  R"doc(
Add temperature model from Stokes coefficient data.

This method adds a temperature model to the coma atmosphere settings using precomputed
Stokes (spherical harmonic) coefficient data. The temperature model uses interpolation
of the Stokes coefficients to compute temperature as a function of position and time.

.. note:: The temperature data type (polynomial or Stokes) must match the density data type.

Parameters
----------
stokes_data : ComaStokesDataset
    Stokes coefficient dataset for temperature distribution.

max_degree : int, default = -1
    Maximum spherical harmonic degree for temperature calculations. Set to -1 to use
    the maximum degree available in the dataset.

max_order : int, default = -1
    Maximum spherical harmonic order for temperature calculations. Set to -1 to use
    the maximum order available in the dataset.

gamma : float, default = 1.33
    Heat capacity ratio (gamma = Cp/Cv) for the gas species. Default value 1.33 is
    appropriate for water vapor.

Examples
--------
.. code-block:: python

  # Create coma settings with Stokes density data
  coma_settings = environment_setup.atmosphere.coma_model_from_stokes_data(
      stokes_data=density_stokes_data,
      molecular_weight=0.018015)

  # Add temperature model with Stokes data
  coma_settings.add_temperature_model(
      stokes_data=temperature_stokes_data,
      max_degree=10,
      max_order=10,
      gamma=1.33)

)doc" );

    m.def( "mars_dtm",
       &tss::marsDtmAtmosphereSettings,
       R"doc(

Function for creating Mars DTM atmospheric settings.

Creates settings for the Mars DTM semiempirical thermosphere model, which is based on the DTM
algorithm originally developed for Earth's thermosphere and adapted for Mars by Bruinsma and
Lemoine (2002). The model reproduces observed densities with approximately 35% uncertainty
(1-σ) outside dust storm periods, with uncertainty increasing by roughly a factor of two
during dust storms.

Bruinsma, S., and F. G. Lemoine (2002), "A preliminary semiempirical thermosphere model of
Mars: DTM-Mars", *Journal of Geophysical Research*, 107(E10), 5085,
doi:10.1029/2001JE001508.

Parameters
----------
space_weather_file : str, default=""
    Path to file containing space weather data for Mars. If an empty string is provided
    (default), the model uses the standard coefficients from Bruinsma & Lemoine (2002)
    without additional space weather corrections. Users can provide a custom file path
    to use updated or mission-specific space weather data.

Returns
-------
AtmosphereSettings
    Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` class

Examples
--------
In this example, we create Mars DTM atmosphere settings with default coefficients:

.. code-block:: python

   # Create Mars DTM atmosphere using default Bruinsma & Lemoine (2002) coefficients
   body_settings.get("Mars").atmosphere_settings = environment_setup.atmosphere.mars_dtm()

In this example, we create Mars DTM atmosphere settings with a custom space weather file:

.. code-block:: python

   # Create Mars DTM atmosphere with custom space weather data
   body_settings.get("Mars").atmosphere_settings = environment_setup.atmosphere.mars_dtm(
       space_weather_file="path/to/custom_mars_space_weather.dat"
   )


)doc" );

#if TUDAT_BUILD_WITH_MCD
    // Factory function for MCD atmosphere
    m.def( "mcd",
           &tss::mcdAtmosphereSettings,
           py::arg( "mcd_data_path" ) = "",
           py::arg( "dust_scenario" ) = 1,
           py::arg( "perturbation_key" ) = 0,
           py::arg( "perturbation_seed" ) = 0.0,
           py::arg( "gravity_wave_length" ) = 0.0,
           py::arg( "high_resolution_mode" ) = 0,
           R"doc(

 Function for creating Mars Climate Database atmosphere model settings.

 Function for settings object, defining atmosphere model using the Mars Climate Database (MCD).
 The MCD provides realistic atmospheric conditions for Mars based on GCM simulations.

 .. note:: **Altitude Convention**: The MCD model expects altitude as "height above areoid" (MOLA zero datum).
           Tudat's flight conditions compute altitude as "height above the shape model surface", which for Mars
           is typically an oblate spheroid that closely approximates the areoid. The difference between these
           definitions is typically less than 0.1%, which is negligible compared to atmospheric variability.

           When ``high_resolution_mode=1``, MCD uses high-resolution MOLA topography to provide more accurate
           atmospheric properties at locations with significant terrain variations.


 Parameters
 ----------
 mcd_data_path : str, default = ""
     Path to MCD data files directory (must end with '/'). If empty string (default), 
     uses the compile-time default path set during CMake configuration. The directory 
     should contain subdirectories like 'clim_aveEUV/', 'dust_high_resol/', etc.

 dust_scenario : int, default = 1
     Dust and solar EUV scenario (controls atmospheric opacity and solar forcing):
     
     **Climatology scenarios (typical Mars year):**
     
     - 1: Climatology, average solar EUV
     - 2: Climatology, minimum solar EUV  
     - 3: Climatology, maximum solar EUV
     
     **Dust storm scenarios (constant opacity=4):**
     
     - 4: Dust storm, minimum solar EUV
     - 5: Dust storm, average solar EUV
     - 6: Dust storm, maximum solar EUV
     
     **Extreme scenarios:**
     
     - 7: Warm scenario (dustier than MY24, maximum solar EUV)
     - 8: Cold scenario (clearer than MY24, minimum solar EUV)
     
     **Mars Year-specific scenarios (with associated solar EUV):**
     
     - 24-35: Mars Years 24 through 35

 perturbation_key : int, default = 0
     Type of atmospheric perturbations to apply:
     
     - 0: None (mean atmospheric state only)
     - 1: Large scale perturbations only (requires ``perturbation_seed``)
     - 2: Small scale (gravity wave) perturbations only (requires ``perturbation_seed`` and ``gravity_wave_length``)
     - 3: Both large and small scale perturbations (requires ``perturbation_seed`` and ``gravity_wave_length``)
     - 4: Both large and small scale perturbations (requires ``perturbation_seed`` and ``gravity_wave_length``)
     - 5: Add n-sigma perturbations, where ``perturbation_seed`` is the multiplier for standard deviations
     
     .. warning:: Perturbations introduce stochastic variability. For reproducible results, use the same
                  ``perturbation_seed`` value. Changing the seed between calls triggers regeneration of
                  perturbation fields.

 perturbation_seed : float, default = 0.0
     Random seed or scaling factor for perturbations (interpretation depends on ``perturbation_key``):
     
     - For ``perturbation_key`` in {1, 2, 3, 4}: Random seed for stochastic perturbations. 
       Changing this value between calls triggers regeneration of perturbation fields.
     - For ``perturbation_key=5``: Multiplier for standard deviations (must be in [-4, 4]).
       Atmospheric properties will be mean ± (``perturbation_seed`` × standard_deviation).

 gravity_wave_length : float, default = 0.0
     Vertical wavelength λ of gravity wave perturbations in meters (required if ``perturbation_key`` 
     in {2, 3, 4}). The horizontal wavelength is automatically set to 10×λ. If set to 0.0, MCD uses 
     its default value of 16000 m. Typical range: 5000-30000 m.

 high_resolution_mode : int, default = 0
     Flag to enable high-resolution topography from MOLA data:
     
     - 0: Use GCM grid resolution (~5.6° × 5.6°)
     - 1: Use high-resolution MOLA topography for improved accuracy near terrain features
     
     .. note:: High-resolution mode is recommended for low-altitude flight near mountains, craters, 
               or other significant terrain features (e.g., Olympus Mons, Valles Marineris).

 Returns
 -------
 AtmosphereSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.atmosphere.AtmosphereSettings` class


 Examples
 --------
 **Example 1**: Basic MCD atmosphere with default settings (climatology, no perturbations):

 .. code-block:: python

    # Create default MCD atmosphere for Mars
    body_settings.get("Mars").atmosphere_settings = environment_setup.atmosphere.mcd()

 **Example 2**: Dust storm scenario with high-resolution topography:

 .. code-block:: python

    # Create MCD atmosphere with dust storm conditions
    body_settings.get("Mars").atmosphere_settings = environment_setup.atmosphere.mcd(
        dust_scenario=5,           # Dust storm with average solar EUV
        high_resolution_mode=1     # Enable high-res MOLA topography
    )

 **Example 3**: Mars Year 32 with stochastic perturbations:

 .. code-block:: python

    # Create MCD atmosphere for Mars Year 32 with gravity wave perturbations
    body_settings.get("Mars").atmosphere_settings = environment_setup.atmosphere.mcd(
        dust_scenario=32,           # Mars Year 32 data
        perturbation_key=3,         # Large + small scale perturbations
        perturbation_seed=42.0,     # Fixed seed for reproducibility
        gravity_wave_length=16000.0 # 16 km vertical wavelength
    )

 **Example 4**: Custom data path and n-sigma perturbations:

 .. code-block:: python

    # Use custom MCD data location and add 2-sigma perturbations
    body_settings.get("Mars").atmosphere_settings = environment_setup.atmosphere.mcd(
        mcd_data_path="/path/to/mcd/data/",
        dust_scenario=1,
        perturbation_key=5,         # n-sigma perturbations
        perturbation_seed=2.0       # Add 2× standard deviation
    )


 See Also
 --------
 :func:`~tudatpy.dynamics.environment_setup.atmosphere.mars_dtm` : Alternative Mars atmosphere model (thermosphere only)
 :func:`~tudatpy.dynamics.environment_setup.atmosphere.exponential_predefined` : Simple exponential Mars atmosphere


 References
 ----------
 .. [1] Millour, E., et al. (2015). "The Mars Climate Database (MCD version 5.2)."
        European Planetary Science Congress, Vol. 10, EPSC2015-438.
 .. [2] Forget, F., et al. (1999). "Improved general circulation models of the Martian atmosphere 
        from the surface to above 80 km." Journal of Geophysical Research, 104(E10), 24155-24175.


     )doc" );
#endif  // TUDAT_BUILD_WITH_MCD
}

}  // namespace atmosphere
}  // namespace environment_setup
}  // namespace dynamics
}  // namespace tudatpy
