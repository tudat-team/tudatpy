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
#include "expose_ground_station.h"

#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
// #include <pybind11/native_enum.h>
#include <tudat/simulation/environment_setup/createGroundStations.h>
#include "tudat/simulation/environment_setup/defaultGroundStationSettings.h"
#include <tudat/simulation/environment_setup/defaultBodies.h>
#include "tudat/simulation/environment_setup/createStationReferencePoint.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace tcc = tudat::coordinate_conversions;
namespace tba = tudat::basic_astrodynamics;

namespace tudatpy
{
namespace dynamics
{
namespace environment_setup
{
namespace ground_station
{

//! Register ground-station construction functions and station-motion settings.
void expose_ground_station_setup( py::module& m )
{
    // Ground station motion settings
    py::enum_< tss::StationMotionModelTypes >( m, "StationMotionModelTypes", R"doc(

         Enumeration of ground station displacement models.

      )doc" )
            .value( "linear", tss::StationMotionModelTypes::linear_station_motion )
            .value( "piecewise_constant", tss::StationMotionModelTypes::piecewise_constant_station_motion )
            .value( "custom", tss::StationMotionModelTypes::custom_station_motion )
            .value( "station_history", tss::StationMotionModelTypes::station_history_motion )
            .value( "body_deformation", tss::StationMotionModelTypes::body_deformation_station_motion )
            .value( "bodycentric_to_barycentric_station_position_motion",
                    tss::StationMotionModelTypes::bodycentric_to_barycentric_station_position_motion );

    py::enum_< tudat::ground_stations::RelativisticStationCoordinateCompatibility >( m, "RelativisticStationCoordinateCompatibility", R"doc(

Enumeration of the coordinate scalings of the relativistic station position correction (see
:func:`relativistic_station_motion`).

      )doc" )
            .value( "unscaled_tcg_to_tcb",
                    tudat::ground_stations::unscaled_tcg_to_tcb_station_coordinates,
                    "Transformation between unscaled (TCG- and TCB-compatible) coordinates, Eq. (11.18) of the IERS Conventions 2010." )
            .value( "tt_to_tdb",
                    tudat::ground_stations::tt_to_tdb_station_coordinates,
                    "Transformation from TT-compatible station coordinates to TDB-compatible coordinates, Eq. (11.19) of the IERS "
                    "Conventions 2010." );

    py::class_< tss::GroundStationMotionSettings, std::shared_ptr< tss::GroundStationMotionSettings > >( m,
                                                                                                         "GroundStationMotionSettings",
                                                                                                         R"doc(

         Base class for providing settings for the motion of a single ground station.

         Non-functional base class for settings for the motion of a single ground station
         Station motion settings requiring additional information must be defined using an object derived from this class.





      )doc" )
            .def_property_readonly( "model_type", &tss::GroundStationMotionSettings::getModelType, R"doc(

         **read-only**

         Type of ground station motion model defined by these settings.

         :type: StationMotionModelTypes

      )doc" );

    m.def( "relativistic_station_motion",
           &tss::relativisticStationMotionSettings,
           py::arg( "external_bodies" ) = std::vector< std::string >{ "Sun" },
           py::arg( "coordinate_compatibility" ) = tudat::ground_stations::tt_to_tdb_station_coordinates,
           R"doc(

 Function for creating settings for the relativistic correction of a ground station position.

 Function for creating settings for the transformation of the station position from the geocentric to the barycentric
 reference system, using the spatial part of Eqs. (11.18) and (11.19) of the IERS Conventions 2010, evaluated at the
 TDB epoch of the host body's centre. The correction is applied to the station position including all other station
 motion and deformation models. The terms due to the simultaneity of events and to the acceleration of the host body
 are omitted. The velocity of the station is the derivative of this correction with respect to TDB.

 A ground station with this motion model must be created with the other bodies of the environment (which is the case
 when using :func:`~tudatpy.dynamics.environment_setup.create_system_of_bodies`, or the ``bodies`` argument of
 :func:`~tudatpy.dynamics.environment_setup.add_ground_station`). The global frame origin must be the solar system
 barycenter.


 Parameters
 ----------
 external_bodies : list[str], default = ["Sun"]
     Bodies (other than the station's host body) whose gravitational potential is included in the correction. An
     empty list gives no external potential.
 coordinate_compatibility : RelativisticStationCoordinateCompatibility, default = RelativisticStationCoordinateCompatibility.tt_to_tdb
     Coordinate scaling of the correction.

 Returns
 -------
 GroundStationMotionSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.ground_station.GroundStationMotionSettings` class
     defining the relativistic station position correction.

 )doc" );

    py::class_< tss::LinearGroundStationMotionSettings,
                std::shared_ptr< tss::LinearGroundStationMotionSettings >,
                tss::GroundStationMotionSettings >( m,
                                                    "LinearGroundStationMotionSettings",
                                                    R"doc(

         Class for defining linear motion (in an Earth-fixed frame) in time of a ground station.

         `GroundStationMotionSettings` derived class for time-linear station motion




      )doc" )
            .def( py::init< const Eigen::Vector3d, const double >( ),
                  py::arg( "linear_velocity" ),
                  py::arg( "reference_epoch" ) = 0.,
                  R"doc(Define linear motion settings for ground station

Parameters
----------

linear_velocity:
    Constant velocity of the station in body-fixed reference frame
reference_epoch:
    Epoch at which the position of the station is known
                  )doc" )
            .def_readwrite( "linear_velocity", &tss::LinearGroundStationMotionSettings::linearVelocity_, R"doc(
                Constant velocity of the station in body-fixed reference frame.

                :type: numpy.ndarray([3,1])
                )doc" )
            .def_readwrite( "reference_epoch", &tss::LinearGroundStationMotionSettings::referenceEpoch_, R"doc(
                Epoch at which the position of the station is defined.

                :type: float
                )doc" );

    py::class_< tss::PiecewiseConstantGroundStationMotionSettings,
                std::shared_ptr< tss::PiecewiseConstantGroundStationMotionSettings >,
                tss::GroundStationMotionSettings >( m,
                                                    "PiecewiseConstantGroundStationMotionSettings",
                                                    R"doc(

         Class for defining piecewise-constant position (e.g. instantaneous change in position at given epochs) of a ground station.

         `GroundStationMotionSettings` derived class for piecewise-constant position of a ground station




      )doc" );

    py::class_< tss::BodyDeformationStationMotionSettings,
                std::shared_ptr< tss::BodyDeformationStationMotionSettings >,
                tss::GroundStationMotionSettings >( m,
                                                    "BodyDeformationStationMotionSettings",
                                                    R"doc(
                    Define station motion settings based on body deformation
                    )doc" );

    py::class_< tss::CustomGroundStationMotionSettings,
                std::shared_ptr< tss::CustomGroundStationMotionSettings >,
                tss::GroundStationMotionSettings >( m,
                                                    "CustomGroundStationMotionSettings",
                                                    R"doc(

         Class for defining custom time-dependent motion of a ground station.

         `CustomGroundStationMotionSettings` derived class for custom time-dependent motion of a ground station




      )doc" );

    py::enum_< tss::WeatherDataSettingsType >( m, "WeatherDataSettingsType" )
            .value( "dsn", tss::WeatherDataSettingsType::dsn_weather_data )
            .value( "estrack", tss::WeatherDataSettingsType::estrack_weather_data );

    py::class_< tss::WeatherDataSettings, std::shared_ptr< tss::WeatherDataSettings > >( m,
                                                                                         "WeatherDataSettings",
                                                                                         R"doc(

         Base class for defining tabulated weather data used by a ground station.

      )doc" )
            .def_property_readonly( "weather_data_type", &tss::WeatherDataSettings::getWeatherDataType )
            .def_property_readonly( "interpolator_settings", &tss::WeatherDataSettings::getInterpolatorSettings );

    py::class_< tss::DsnWeatherDataSettings, std::shared_ptr< tss::DsnWeatherDataSettings >, tss::WeatherDataSettings >(
            m,
            "DsnWeatherDataSettings",
            R"doc(

         Settings for DSN tabulated weather data.

      )doc" );

    py::class_< tss::EstrackWeatherDataSettings, std::shared_ptr< tss::EstrackWeatherDataSettings >, tss::WeatherDataSettings >(
            m,
            "EstrackWeatherDataSettings",
            R"doc(

         Settings for ESTRACK tabulated weather data.

      )doc" );

    py::class_< tss::GroundStationSettings, std::shared_ptr< tss::GroundStationSettings > >( m,
                                                                                             "GroundStationSettings",
                                                                                             R"doc(

         Base class for providing settings for the creation of a ground station.


      )doc" )
            .def_property( "station_position",
                           &tss::GroundStationSettings::getGroundStationPosition,
                           &tss::GroundStationSettings::resetGroundStationPosition,
                           R"doc(
                           Position of the ground station in body-fixed frame. The position is interpreted based on the value of the ``position_element_type`` property.

                           :type: numpy.ndarray([3,1])

                           )doc" )

            .def_property_readonly( "station_name", &tss::GroundStationSettings::getStationName, R"doc(
                Name of the ground station.

                :type: str

                )doc" )
            .def_property_readonly( "position_element_type", &tss::GroundStationSettings::getPositionElementType, R"doc(

            Element type of the ground station position, defining the interpretation of the ``station_position`` property.

            :type: ~tudatpy.astro.element_conversion.PositionElementTypes

                )doc" )
            .def_property( "station_motion_settings",
                           &tss::GroundStationSettings::getStationMotionSettings,
                           &tss::GroundStationSettings::setStationMotionSettings,
                           R"doc(

                List of motion settings for the ground station, defining time-variations of the station position.

                :type: list[ GroundStationMotionSettings ]

                    )doc" )
            .def_property( "weather_data_settings",
                           &tss::GroundStationSettings::getWeatherDataSettings,
                           &tss::GroundStationSettings::setWeatherDataSettings,
                           R"doc(

                Settings for tabulated weather data to load when the ground station is created.

                :type: WeatherDataSettings

                    )doc" )
            .def( "add_station_motion_settings",
                  &tss::GroundStationSettings::addStationMotionSettings,
                  py::arg( "station_motion_settings" ),
                  R"doc(Add a motion model to the existing motion models of the station.

 Parameters
 ----------

 station_motion_settings: GroundStationMotionSettings
     GroundStationMotionSettings object defining the motion model to be added

)doc" );

    m.def( "add_motion_model_to_each_ground_station",
           &tss::addStationMotionModelToEachGroundStation,
           py::arg( "ground_station_settings_list" ),
           py::arg( "station_motion_setting" ),
           R"doc(

         Add a station motion model to each ground station in a settings list.

         The motion model is appended to the models already present in each settings object. The supplied ground station
         settings are modified in place.

         Parameters
         ----------
         ground_station_settings_list : list[GroundStationSettings]
             Settings of the ground stations to which the motion model is added.
         station_motion_setting : GroundStationMotionSettings
             Settings defining the motion model to add to each station.

      )doc" );

    m.def( "set_dsn_weather_data_in_ground_station_settings",
           py::overload_cast< const std::vector< std::shared_ptr< tss::GroundStationSettings > >&,
                              const std::vector< std::string >&,
                              std::shared_ptr< tudat::interpolators::InterpolatorSettings >,
                              const std::map< int, std::vector< std::string > >& >( &tss::setDsnWeatherDataInGroundStationSettings ),
           py::arg( "ground_station_settings_list" ),
           py::arg( "weather_file_names" ),
           py::arg_v( "interpolator_settings", tudat::interpolators::linearInterpolation( ), "..." ),
           py::arg_v( "ground_stations_per_complex",
                      tss::getDefaultDsnStationNamesPerComplex( ),
                      "tudatpy.dynamics.environment_setup.ground_station.get_default_dsn_station_names_per_complex()" ),
           R"doc(

 Add DSN weather data settings to a list of ground station settings.

      )doc" );

    m.def( "set_estrack_weather_data_in_ground_station_settings",
           py::overload_cast< const std::vector< std::shared_ptr< tss::GroundStationSettings > >&,
                              const std::vector< std::string >&,
                              const std::string,
                              std::shared_ptr< tudat::interpolators::InterpolatorSettings > >(
                   &tss::setEstrackWeatherDataInGroundStationSettings ),
           py::arg( "ground_station_settings_list" ),
           py::arg( "weather_file_names" ),
           py::arg( "ground_station_name" ),
           py::arg_v( "interpolator_settings",
                      tudat::interpolators::cubicSplineInterpolation(
                              tudat::interpolators::AvailableLookupScheme::huntingAlgorithm,
                              tudat::interpolators::BoundaryInterpolationType::use_boundary_value_with_warning ),
                      "..." ),
           R"doc(

 Add ESTRACK weather data settings to a ground station in a list of ground station settings.

      )doc" );

    m.def( "basic_station",
           &tss::groundStationSettings,
           py::arg( "station_name" ),
           py::arg( "station_nominal_position" ),
           py::arg( "station_position_element_type" ) = tcc::cartesian_position,
           py::arg( "station_motion_settings" ) = std::vector< std::shared_ptr< tss::GroundStationMotionSettings > >( ),
           py::arg( "weather_data_settings" ) = nullptr,
           R"doc(

 Function for creating settings for a ground station

 Function for creating settings for a ground station, defining only its name, body-fixed position, and (optionally) time-variations of its position


 Parameters
 ----------
 station_name : string
     Name (unique identifier) by which the station is to be known.
 station_position_element_type : PositionElementTypes, default = cartesian_position
     Type of elements for ``station_nominal_position``.
 station_nominal_position : numpy.ndarray([3,1])
     Nominal position of the station in a body-fixed frame. Depending on the choice of ``station_position_element_type`` input, this vector must contain

     * Cartesian (for ``cartesian_position`` input): :math:`[x,y,z]`, denoting :math:`x-`, :math:`y-` and :math:`z-` components of body-fixed position (w.r.t body-fixed frame origin, typically center of mass)
     * Spherical (for ``spherical_position`` input): :math:`[r,\phi',\theta]`, denoting distance from body-fixed frame origin (typically center of mass), latitude and longitude
     * Geodetic (for ``geodetic_position`` input): :math:`[h,\phi,\theta]`, denoting the altitude w.r.t. the body shape model, geodetic latitude and longitude. Note that, in this case, the conversion to Cartesian position depends on the body's shape model

 station_motion_settings : list[ GroundStationMotionSettings ], default = None
     List of settings defining time-variations of the individual ground station
 Returns
 -------
 GroundStationSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.ground_station.GroundStationSettings` defining settings of the to be created ground station





 Examples
 --------
 In this example, we create a station using geodetic coordinates at the approximate location of the city of Delft, and no motion settings:

 .. code-block:: python

   # Define the position of the ground station on Earth
   station_altitude = 0.0
   delft_latitude = np.deg2rad(52.00667)
   delft_longitude = np.deg2rad(4.35556)

   # Create ground station settings
   ground_station_settings = environment_setup.ground_station.basic_station(
       "TrackingStation",
        [station_altitude, delft_latitude, delft_longitude],
        element_conversion.geodetic_position_type)

   # Append station settings to existing (default is empty) list
           body_settings.get( "Earth" ).ground_station_settings.append( ground_station_settings )


     )doc" );

    m.def( "dsn_station",
           &tss::getDsnStationSetting,
           py::arg( "station_name" ),
           R"doc(

 Function for creating settings for single DSN station (see :func:`~dsn_stations`)

 Parameters
 ----------
 station_name : str
     String with the name of the station, e.g. "DSS-12"

 Returns
 -------
 GroundStationSettings
     Settings to create DSN station

 Examples
 --------

 In this example, settings for the 70m DSN dishes DSS-14, DSS-43 and DSS-63 are created and assigned to the Earth ground station settings:

 .. code-block:: python

     from tudatpy.dynamics.environment_setup.ground_station import \
         dsn_station

     dss_station_names = ["DSS-14", "DSS-43", "DSS-63"]
     dss_station_settings = [dsn_station(station_name) for station_name in dss_station_names]

     body_settings.get("Earth").ground_station_settings = dss_station_settings


     )doc" );
    m.def( "dsn_stations",
           &tss::getDsnStationSettings,
           R"doc(

 Function for creating settings for all DSN stations

 Function for creating settings for all DSN stations, defined by nominal positions and linear velocities, as defined
 by Cartesian elements in *DSN No. 810-005, 301, Rev. O*,  see `this link <https://deepspace.jpl.nasa.gov/dsndocs/810-005/301/301O.pdf>`__.
 DSS-65 is defined with a piecewise-constant motion, according to `this link <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/stations/a_old_versions/earthstns_itrf93_050714.cmt>`__, as it was moved to a new location in 2005.
 Note that calling these settings will use the Cartesian elements provided in this document (in ITRF93) and apply them to the Earth-fixed
 station positions, regardless of the selected Earth rotation model.

 Returns
 -------
 list[ GroundStationSettings ]
     List of settings to create DSN stations

     )doc" );

    m.def( "evn_stations",
           &tss::getEvnStationSettings,
           R"doc(

 Function for creating settings for all EVN stations.

 Function for creating settings for all EVN stations. EVN stations are defined by nominal positions and linear velocities, as defined by the glo.sit station file, see `this link <https://gitlab.com/gofrito/pysctrack/-/blob/master/cats/glo.sit?ref_type=heads>`__.
 Note that calling these settings will use the Cartesian elements provided by these documents and apply them to the Earth-fixed station positions, regardless of the selected Earth rotation model.

 Returns
 -------
 list[ GroundStationSettings ]
     List of settings to create EVN stations

     )doc" );

    m.def( "radio_telescope_stations",
           &tss::getRadioTelescopeStationSettings,
           R"doc(

 Function for creating settings for all DSN and EVN stations, combining :func:`~evn_stations` and :func:`dsn_stations`

 Returns
 -------
 list[ GroundStationSettings ]
     List of settings to create DSN + EVN stations

     )doc" );

    m.def( "optical_telescope_stations",
           &tss::getMPCStationSettings,
           R"doc(

 Function for creating settings for existing MPC stations.
 Observatories cartesian coordinates were converted using parallax info (rho_cos_phi, rho_sin_phi)
 from Project Pluto's list: https://www.projectpluto.com/mpc_stat.txt,
 which has them in geodetic coordinates. Some station positions (e.g. amateur observatories positions) might lack accuracy,
 especially in altitude, as it is sometimes not clear whether their altitude has to be intended as a geodetic or ellipsoidal altitude.
 For more info and insights, please check the following discussion: https://www.projectpluto.com/mpc_stat.htm.

 Returns
 -------
 list[ GroundStationSettings ]
     List of settings to create MPC stations

     )doc" );

    m.def( "default_ilrs_sinex_state_file",
           &tss::getDefaultIlrsSinexStateFilePath,
           R"doc(

 Return the default SINEX station-state file path for ILRS station loading.

 Returns
 -------
 str
     Path to the default ILRS SINEX state file in Tudat resources.

     )doc" );

    m.def( "default_ilrs_sinex_eccentricity_file",
           &tss::getDefaultIlrsSinexEccentricityFilePath,
           R"doc(

 Return the default SINEX station-eccentricity file path for ILRS station loading.

 Returns
 -------
 str
     Path to the default ILRS SINEX eccentricity file in Tudat resources.

     )doc" );

    m.def( "ilrs_stations_from_sinex_domes",
           &tss::getIlrsStationSettingsFromSinexDomes,
           py::arg( "domes_ids" ),
           py::arg( "sinex_state_file" ),
           py::arg( "sinex_eccentricity_file" ),
           py::arg( "throw_exception_on_missing_data" ) = false,
           R"doc(

 Function for creating ILRS station settings from SINEX state/eccentricity files using IERS DOMES identifiers.

 Parameters
 ----------
 domes_ids : List[str]
     List of IERS DOMES identifiers to load (for example, ``["14201S018"]`` for Wettzell).
 sinex_state_file : str
     Path to SINEX station-state file (position/velocity).
 sinex_eccentricity_file : str
     Path to SINEX eccentricity file. Set to an empty string to disable eccentricity offsets.
 throw_exception_on_missing_data : bool, default = False
     Whether to throw when requested stations cannot be mapped or found.

 Returns
 -------
 list[ GroundStationSettings ]
     Ground station settings for the selected ILRS DOMES stations.

     )doc" );

    m.def( "linear_station_motion",
           &tss::linearGroundStationMotionSettings,
           py::arg( "linear_velocity" ),
           py::arg( "reference_epoch" ) = 0.0,
           R"doc(

 Function for creating settings for a linear station motion

 Function for creating settings for a linear station motion, implementing :math:`\Delta \mathbf{r}=\dot{\mathbf{r}}(t-t_{0})`.


 Parameters
 ----------
 linear_velocity : numpy.ndarray([3,1])
     Linear velocity :math:`\dot{\mathbf{r}}` of the station (in m/s)
 reference_epoch : float, default = 0.0
     Reference epoch :math:`t_{0}`, in seconds since J2000 TDB.
 Returns
 -------
 GroundStationMotionSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.ground_station.GroundStationMotionSettings` derived :class:`~tudatpy.dynamics.environment_setup.ground_station.LinearGroundStationMotionSettings` class



     )doc" );

    m.def( "piecewise_constant_station_motion",
           &tss::piecewiseConstantGroundStationMotionSettings,
           py::arg( "displacement_list" ),
           R"doc(

 Function for creating settings for a piecewise constant ground station position variation

 Function for creating settings for a piecewise constant ground station position. Using this model, the added station velocity in a body-fixed frame :math:`\dot{\mathbf{r}}` is
 always zero, but its displacement :math:`\Delta\mathbf{r}` is set according to the input list, which contains a list of times and displacements :math:`[t_{i},\Delta\mathbf{r}_{i}]`.
 When the resulting model is queried at a given time :math:`t`, the nearest lower neighbour :math:`t_{i}` from this list is found, and the associated :math:`\Delta\mathbf{r}_{i}` is applied.


 Parameters
 ----------
 displacement_list : dict[float, numpy.ndarray([3, 1])]
     Dictionary with epochs :math:`t_{i}` in seconds since J2000 TDB as keys, and the associated displacements :math:`\Delta\mathbf{r}_{i}` as values.
 Returns
 -------
 GroundStationMotionSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.ground_station.GroundStationMotionSettings` derived :class:`~tudatpy.dynamics.environment_setup.ground_station.PiecewiseConstantGroundStationMotionSettings` class






     )doc" );

    m.def( "body_deformation_station_motion",
           &tss::bodyDeformationStationMotionSettings,
           py::arg( "throw_exception_when_not_available" ) = true,
           R"doc(
Apply the body's configured shape-deformation models to this station or surface reference point.

The displacement is evaluated at the station epoch in TDB and added in body-fixed axes. Configure the physical models
through the body's shape_deformation_settings; this function connects those existing models to the station.
Keep throw_exception_when_not_available true when every configured contribution is required by the analysis.
           )doc" );

    m.def( "custom_station_motion",
           &tss::customGroundStationMotionSettings,
           py::arg( "custom_displacement_function" ),
           R"doc(

 Function for creating settings for a custom ground station position variation

 Function for creating settings for a custom ground station position. An arbitrary user-defined function of the signature :math:`\Delta\mathbf{r}=\Delta\mathbf{r}(t)` is provided and
 applied to the station position


 Parameters
 ----------
 custom_displacement_function : callable[[float], numpy.ndarray([3, 1])]
     Function returning :math:`\Delta\mathbf{r}`, with time :math:`t` in seconds since J2000 TDB as input.
 Returns
 -------
 GroundStationMotionSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.ground_station.GroundStationMotionSettings` derived :class:`~tudatpy.dynamics.environment_setup.ground_station.CustomGroundStationMotionSettings` class






     )doc" );

    m.def( "get_approximate_dsn_ground_station_positions", &tss::getApproximateDsnGroundStationPositions, R"doc(

    This function returns the approximate positions of the DSN ground stations.

    The function returns ground station positions for all present ground stations specified in `DSN 810-005, 301 Coverage and Geometry, Revision O (2024), DSN/JPL <https://deepspace.jpl.nasa.gov/dsndocs/810-005/301/301O.pdf>`__. Additionally, historic positions for DSS-12, DSS-42 and DSS-61 are provided, retrieved from `NAIF <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/stations/a_old_versions/dsnstns.cmt>`__. The positions of the ground stations are specified at 2003.0 with respect to ITRF93.

    Returns
    -------
    dict[str, numpy.ndarray([3,1])]
        Dictionary mapping DSN station names (str, format: "DSS-<id>") to approximate positions.

        )doc" );

    m.def( "get_default_dsn_station_names_per_complex", &tss::getDefaultDsnStationNamesPerComplex, R"doc(

    Return the default station names for each DSN complex.

    Returns
    -------
    dict[int, list[str]]
        Dictionary mapping each DSN complex identifier to its station names.

        )doc" );

    m.def( "get_radio_telescope_positions", &tss::getCombinedApproximateGroundStationPositions, R"doc(

    This function returns the positions of DSN ground stations and VLBI stations.

    The function returns ground station positions for all DSN ground stations (as given by :func:`~tudatpy.dynamics.environment_setup.ground_station.get_approximate_dsn_ground_station_positions`) and VLBI stations (as given by :func:`~tudatpy.dynamics.environment_setup.ground_station.get_vlbi_station_positions`).

    Returns
    -------
    dict[str, numpy.ndarray([3,1])]
        Dictionary mapping station name to positions.

        )doc" );

    m.def( "approximate_ground_stations_position", &tss::getCombinedApproximateGroundStationPositions, R"doc(

         Retrieve approximate body-fixed positions of DSN and VLBI stations.

         This function returns the same station catalogue as :func:`~get_radio_telescope_positions`, combining the positions
         returned by :func:`~get_approximate_dsn_ground_station_positions` and :func:`~get_vlbi_station_positions`.

         Returns
         -------
         dict[str, numpy.ndarray[numpy.float64[3, 1]]]
             Dictionary mapping each station name to its Earth-fixed Cartesian position, in m.

      )doc" );

    m.def( "get_vlbi_station_positions", &tss::getVlbiStationPositions, R"doc(

    This function returns the positions of VLBI stations.

    The VLBI station positions are retrieved from `pysctrack <https://gitlab.com/gofrito/pysctrack/-/raw/master/cats/glo.sit>`__.

    Returns
    -------
    dict[str, numpy.ndarray([3,1])]
        Dictionary mapping station name to positions.

        )doc" );
    m.def( "get_vlbi_station_velocities", &tss::getVlbiStationVelocities, R"doc(

    This function returns the velocities of VLBI stations.

    The VLBI station velocities are retrieved from `pysctrack <https://gitlab.com/gofrito/pysctrack/-/raw/master/cats/glo.sit>`__.

    Returns
    -------
    dict[str, numpy.ndarray([3,1])]
        Dictionary mapping station name to velocities.

        )doc" );

    // Build source-profile trajectories from neutral input products, then bind absolute station motion.
    m.def( "ign_station_history",
           &tudat::ground_stations::createIgnStationCoordinateHistory,
           py::arg( "solutions" ),
           py::arg( "discontinuities" ),
           py::arg( "site" ),
           py::arg( "post_seismic" ) = nullptr,
           py::arg( "seasonal" ) = nullptr,
           py::arg( "enu_to_terrestrial" ) = py::none( ),
           py::arg( "domes_selector" ) = py::none( ),
           py::arg( "velocity_consistency_tolerance_metres_per_julian_year" ) = 1.0e-8,
           R"doc(

 Function for creating a station history from IGN (ITRF) station products.

 Function for creating the station history of a single station from the SINEX products of an ITRF solution, as
 distributed by IGN: the station coordinate solution, the discontinuity file and, optionally, the post-seismic
 deformation and seasonal-term files. The station position is the sum of:

 * the linear (position and velocity) solution that is valid at the epoch. Each position solution must lie within a
   single velocity interval of the discontinuity file, and solutions sharing a velocity interval must have consistent
   velocities;
 * all post-seismic deformation terms (logarithmic and exponential) of earthquakes before the epoch. An amplitude with
   parameter index :math:`i` is paired with the relaxation time with index :math:`i+1`; relaxation times are
   interpreted in Julian years (although labelled ``m`` in the IGN files);
 * the annual and semi-annual seasonal terms of the valid position solution, with a phase given by the fraction of the
   actual calendar year.

 The time conventions of the IGN products are used: elapsed time for the linear and post-seismic terms is the
 difference of UTC modified Julian dates divided by 365.25, without leap seconds. Solution validity intervals are
 :math:`[t_{start}, t_{end})` in TAI, so a new solution applies from its start epoch; epochs in a gap between
 solutions, or outside them, raise an exception.

 The ``frame``, ``release``, ``origin`` and ``tide_system`` attributes of all products must be set by the user, and be
 identical; they are not read from the files.


 Parameters
 ----------
 solutions : SinexStationHistory
     Station coordinate solution (for instance ``ITRF2020-IGS-TRF.SNX``).
 discontinuities : SinexStationHistory
     Discontinuity file, giving the validity intervals of the position and velocity solutions.
 site : SinexSiteKey
     Site code and point code of the station.
 post_seismic : SinexStationHistory, default = None
     Post-seismic deformation parameters.
 seasonal : SinexStationHistory, default = None
     Seasonal-term parameters.
 enu_to_terrestrial : numpy.ndarray[numpy.float64[3, 3]], default = None
     Constant rotation matrix from the local (east, north, up) basis to body-fixed Cartesian coordinates. Required if
     the post-seismic or seasonal terms are given in east, north, up components.
 domes_selector : str, default = None
     DOMES number of the monument, required when several monuments share the site and point code.
 velocity_consistency_tolerance_metres_per_julian_year : float, default = 1e-8
     Tolerance on the difference between velocity components of solutions that share a velocity interval.

 Returns
 -------
 StationCoordinateHistory
     Station coordinate history of the station marker.


 Examples
 --------
 .. code-block:: python

     from tudatpy.data_input.environment_data import station_history
     from tudatpy.dynamics.environment_setup import ground_station

     history = ground_station.ign_station_history(
         coordinates, discontinuities, station_history.SinexSiteKey( "AB07", "A" ),
         post_seismic=psd, enu_to_terrestrial=enu_basis )
     station_settings = ground_station.basic_station(
         "AB07", nominal_position,
         station_motion_settings=[ ground_station.station_history_motion( history ) ] )

 )doc" );
    m.def( "transform_station_history",
           &tudat::ground_stations::transformStationCoordinateHistory,
           py::arg( "history" ),
           py::arg( "catalogue" ),
           py::arg( "target_frame" ),
           py::arg( "target_origin" ),
           py::arg( "target_tide_system" ),
           py::arg( "acknowledge_aligned_releases" ) = false,
           R"doc(

 Function for transforming a station coordinate history to another terrestrial reference frame realization.

 Function for transforming a station coordinate history to another terrestrial reference frame realization, using a
 table of linearized Helmert transformation parameters. The transformation is applied once to the coefficients, and a
 new history is returned: each linear solution is transformed at its reference epoch (position, and velocity including
 the rates of the transformation parameters), and the post-seismic and seasonal amplitudes are rotated and scaled.
 Compared with evaluating the transformation at each epoch, the position differs by up to a few micrometres for the
 oldest ITRF realizations (because the rates of the parameters are per calendar year, and those of the station per
 Julian year). The reference origin and tide system are not converted, and must be equal to those of the input history.


 Parameters
 ----------
 history : StationCoordinateHistory
     Station coordinate history that is to be transformed.
 catalogue : TerrestrialFrameTransformationData
     Table of transformation parameters (see :ref:`terrestrial_frame_data`).
 target_frame : TerrestrialFrameIdentity
     Realization and release of the target frame.
 target_origin : str
     Reference origin of the output (must be equal to that of the input).
 target_tide_system : str
     Tide system of the output (must be equal to that of the input).
 acknowledge_aligned_releases : bool, default = False
     Whether a release that differs from that of the catalogue (for instance an update of a realization) is accepted as
     aligned with it. This does not add a transformation.

 Returns
 -------
 StationCoordinateHistory
     Station coordinate history in the target frame.


 Examples
 --------
 .. code-block:: python

     from tudatpy.data_input.environment_data import terrestrial_frames
     from tudatpy.dynamics.environment_setup import ground_station

     target = terrestrial_frames.TerrestrialFrameIdentity( )
     target.realization = target.release = "ITRF2014"
     transformed = ground_station.transform_station_history(
         history, terrestrial_frames.official_itrf2020_parameters( ), target,
         target_origin=history.metadata.origin, target_tide_system=history.metadata.tide_system )

 )doc" );
    py::class_< tss::StationHistoryGroundStationMotionSettings,
                std::shared_ptr< tss::StationHistoryGroundStationMotionSettings >,
                tss::GroundStationMotionSettings >( m, "StationHistoryGroundStationMotionSettings", R"doc(

 Class for defining settings for station motion from a station history.

 `GroundStationMotionSettings` derived class for station motion from a station history. Instances of this class are
 created by the :func:`station_history_motion` function.

      )doc" );
    m.def( "station_history_motion",
           &tss::stationHistoryGroundStationMotionSettings,
           py::arg( "history" ),
           py::arg( "reference_point_offsets" ) =
                   std::vector< std::shared_ptr< tudat::ground_stations::StationReferencePointOffsetHistory > >( ),
           py::arg( "source_time_function" ) = py::none( ),
           R"doc(

 Function for creating settings for station motion from a station coordinate history.

 Function for creating settings for station motion from a station coordinate history, optionally moved to another point
 of the station by reference point offsets (for instance from the marker to the antenna reference point). The nominal
 position of the ground station is subtracted, so that the resulting station position is equal to the position from the
 history plus the offsets. By default, the epoch is converted to the time coordinates of the station history using the
 IGN conventions (see :func:`~tudatpy.dynamics.environment.ign_nominal_utc_station_time_from_tdb`), with the time scale
 converter of the body's rotation model when available.


 Parameters
 ----------
 history : StationCoordinateHistory
     Station coordinate history (for instance from :func:`ign_station_history`).
 reference_point_offsets : list[StationReferencePointOffsetHistory], default = []
     Offset histories added to the position of the history, in order: the first starts at the point of the history, and
     each next one at the end point of the previous one. They must be given in the frame of the history.
 source_time_function : callable[[astro.time_representation.Time, bool], StationHistoryTime], default = None
     Custom conversion of the TDB epoch (first argument) to the time coordinates of the station history, setting the
     rates when the second argument is true. Its time convention must match that of the station history.

 Returns
 -------
 StationHistoryGroundStationMotionSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.ground_station.GroundStationMotionSettings` derived
     :class:`~tudatpy.dynamics.environment_setup.ground_station.StationHistoryGroundStationMotionSettings` class.

 )doc" );

    // Bind reference-point geometry and ordered offset composition before constructing a named station.
    m.def( "geocentric_reference_point_geometry",
           &tudat::ground_stations::createGeocentricStationReferencePointGeometry,
           py::arg( "anchor_metres" ),
           py::arg( "frame" ),
           R"doc(

 Function for creating a local basis with a geocentric vertical, used to convert local offsets of station reference
 points to Cartesian coordinates. The basis is fixed at the given position.


 Parameters
 ----------
 anchor_metres : numpy.ndarray[numpy.float64[3, 1]]
     Cartesian position at which the basis is defined, in metres.
 frame : str
     Terrestrial reference frame of the position.

 Returns
 -------
 StationReferencePointGeometry
     Local basis, with the up direction along the geocentric radial direction.

 )doc" );
    m.def( "geodetic_reference_point_geometry",
           &tudat::ground_stations::createGeodeticStationReferencePointGeometry,
           py::arg( "anchor_metres" ),
           py::arg( "frame" ),
           py::arg( "equatorial_radius_metres" ),
           py::arg( "flattening" ),
           R"doc(

 Function for creating a local basis with a geodetic vertical, used to convert local offsets of station reference
 points to Cartesian coordinates. The basis is fixed at the given position.


 Parameters
 ----------
 anchor_metres : numpy.ndarray[numpy.float64[3, 1]]
     Cartesian position at which the basis is defined, in metres.
 frame : str
     Terrestrial reference frame of the position.
 equatorial_radius_metres : float
     Equatorial radius of the reference ellipsoid, in metres.
 flattening : float
     Flattening of the reference ellipsoid.

 Returns
 -------
 StationReferencePointGeometry
     Local basis, with the up direction normal to the reference ellipsoid.

 )doc" );
    m.def( "sinex_antenna_reference_point_offsets",
           &tudat::ground_stations::createSinexAntennaReferencePointOffsets,
           py::arg( "marker" ),
           py::arg( "equipment" ),
           py::arg( "reference_point_name" ),
           py::arg( "geometry" ) = py::none( ),
           py::arg( "header_context" ) = py::none( ),
           py::arg( "domes_selector" ) = py::none( ),
           py::arg( "equipment_site" ) = py::none( ),
           R"doc(

 Function for creating the offsets from the marker to the antenna reference point (ARP) of a GNSS station.

 Function for creating the offsets from the marker to the antenna reference point of a GNSS station, from the
 marker-to-ARP eccentricities of a SINEX file (for instance the IGS station log summary). The offsets are added to the
 station coordinate history of the marker with :func:`station_history_motion` or :func:`reference_point_station`.
 Each eccentricity is associated with the antenna that is valid on its interval; the antenna type, radome and serial
 number are retained. Eccentricities in up, north, east components are converted with the given local basis. The phase
 centre is not modelled. At an epoch where the antenna changes, the new eccentricity is used.


 Parameters
 ----------
 marker : StationHistoryMetadata
     Description of the station coordinate history of the marker (for instance ``history.metadata``); the offsets start
     at this marker and are given in its frame.
 equipment : SinexStationHistory
     SINEX file with the antenna and eccentricity histories of the station.
 reference_point_name : str
     Name of the antenna reference point (also used as station name by :func:`reference_point_station`).
 geometry : StationReferencePointGeometry, default = None
     Local basis, required for eccentricities in up, north, east components.
 header_context : SinexEquipmentHeaderTimeContext, default = None
     Validity interval to use for records whose start or end refers to a file header epoch that is not given.
 domes_selector : str, default = None
     DOMES number of the monument, required when several monuments share the site and point code.
 equipment_site : SinexSiteKey, default = None
     Site and point code of the station in the equipment file, if it differs from that of the marker (for instance
     ``abmf`` instead of ``ABMF``). It must have the same DOMES number.

 Returns
 -------
 StationReferencePointOffsetHistory
     Offsets from the marker to the antenna reference point.

 )doc" );
    m.def( "reference_point_station",
           &tss::stationReferencePointSettings,
           py::arg( "history" ),
           py::arg( "reference_point_offsets" ),
           py::arg( "nominal_position" ),
           py::arg( "source_time_function" ) = py::none( ),
           R"doc(

 Function for creating settings for a ground station located at a reference point of a station.

 Function for creating settings for a ground station whose position is given by a station coordinate history plus
 reference point offsets (for instance to the antenna reference point, see :func:`sinex_antenna_reference_point_offsets`).
 The name of the end point of the last offsets is used as station name, and :func:`station_history_motion` is used as
 station motion.


 Parameters
 ----------
 history : StationCoordinateHistory
     Station coordinate history (for instance of the station marker).
 reference_point_offsets : list[StationReferencePointOffsetHistory]
     Offset histories from the point of the history to the reference point, in order.
 nominal_position : numpy.ndarray[numpy.float64[3, 1]]
     Nominal Cartesian position of the station, in the body-fixed frame.
 source_time_function : callable[[astro.time_representation.Time, bool], StationHistoryTime], default = None
     Custom conversion of the epoch to the time coordinates of the station history (see :func:`station_history_motion`).

 Returns
 -------
 GroundStationSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.ground_station.GroundStationSettings` defining the
     ground station.

 )doc" );
}

}  // namespace ground_station
}  // namespace environment_setup
}  // namespace dynamics
}  // namespace tudatpy
