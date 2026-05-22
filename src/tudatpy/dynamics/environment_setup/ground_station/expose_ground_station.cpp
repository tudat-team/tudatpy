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

void expose_ground_station_setup( py::module& m )
{
    // Ground station motion settings
    py::enum_< tss::StationMotionModelTypes >( m, "StationMotionModelTypes" )
            .value( "linear", tss::StationMotionModelTypes::linear_station_motion )
            .value( "piecewise_constant", tss::StationMotionModelTypes::piecewise_constant_station_motion )
            .value( "custom", tss::StationMotionModelTypes::custom_station_motion )
            .value( "body_deformation", tss::StationMotionModelTypes::body_deformation_station_motion )
            .value( "bodycentric_to_barycentric_station_position_motion",
                    tss::StationMotionModelTypes::bodycentric_to_barycentric_station_position_motion );

    py::class_< tss::GroundStationMotionSettings, std::shared_ptr< tss::GroundStationMotionSettings > >( m,
                                                                                                         "GroundStationMotionSettings",
                                                                                                         R"doc(

         Base class for providing settings for the motion of a single ground station.

         Non-functional base class for settings for the motion of a single ground station
         Station motion settings requiring additional information must be defined using an object derived from this class.





      )doc" )
            .def_property_readonly( "model_type", &tss::GroundStationMotionSettings::getModelType );

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

                  :param linear_velocity: Constant velocity of the station in body-fixed reference frame
                  :param reference_epoch: Epoch at which the position of the station is known
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
                    )doc" )
            .def( py::init< const bool >( ), py::arg( "fail_if_not_available" ) = true );

    py::class_< tss::CustomGroundStationMotionSettings,
                std::shared_ptr< tss::CustomGroundStationMotionSettings >,
                tss::GroundStationMotionSettings >( m,
                                                    "CustomGroundStationMotionSettings",
                                                    R"doc(

         Class for defining custom time-dependent motion of a ground station.

         `CustomGroundStationMotionSettings` derived class for custom time-dependent motion of a ground station




      )doc" );

    py::class_< tss::GroundStationSettings, std::shared_ptr< tss::GroundStationSettings > >( m,
                                                                                             "GroundStationSettings",
                                                                                             R"doc(

         Base class for providing settings for the creation of a ground station.


      )doc" )
            .def_property( "station_position",
                           &tss::GroundStationSettings::getGroundStationPosition,
                           &tss::GroundStationSettings::resetGroundStationPosition )

            .def_property_readonly( "station_name", &tss::GroundStationSettings::getStationName )
            .def_property_readonly( "position_element_type", &tss::GroundStationSettings::getPositionElementType )
            .def_property_readonly( "station_motion_settings", &tss::GroundStationSettings::getStationMotionSettings );

    m.def( "add_motion_model_to_each_ground_station",
           &tss::addStationMotionModelToEachGroundStation,
           py::arg( "ground_station_settings_list" ),
           py::arg( "station_motion_setting" ) );

    m.def( "basic_station",
           &tss::groundStationSettings,
           py::arg( "station_name" ),
           py::arg( "station_nominal_position" ),
           py::arg( "station_position_element_type" ) = tcc::cartesian_position,
           py::arg( "station_motion_settings" ) = std::vector< std::shared_ptr< tss::GroundStationMotionSettings > >( ),
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
     * Cartesian (for ``cartesian_position`` input) - :math:`[x,y,z]`, denoting :math:`x-`, :math:`y-` and :math:`z-` components of body-fixed position (w.r.t body-fixed frame origin, typically center of mass)
     * Spherical (for ``spherical_position`` input) - - :math:`[r,\phi',\theta]`, denoting distance from body-fixed frame origin (typically center of mass), latitude and longitude
     * Geodetic (for ``geodetic_position`` input) - - :math:`[h,\phi,\theta]`, denoting the altitude w.r.t. the body shape model, geodetic latitude and longitude. Note that, in this case, the conversion to Cartesian position depends on the body's shape model
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
 reference_epoch : astro.time_representation.Time, default = 0.0
     Reference epoch :math:`t_{0}` (Time object representing seconds since J2000 TDB)
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
 displacement_list : dict[astro.time_representation.Time,numpy.ndarray([3,1])]
     Dictionary with the epochs :math:`t_{i}` as keys (as Time objects), and the associated displacement :math:`\Delta\mathbf{r}_{i}` as values
 Returns
 -------
 GroundStationMotionSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.ground_station.GroundStationMotionSettings` derived :class:`~tudatpy.dynamics.environment_setup.ground_station.PiecewiseConstantGroundStationMotionSettings` class






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
 custom_displacement_function : callable[[:class:`~tudatpy.astro.time_representation.Time`],numpy.ndarray([3,1])]
     Function returning :math:`\Delta\mathbf{r}`, with the time :math:`t` (as Time object) as input.
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

    m.def( "get_radio_telescope_positions", &tss::getCombinedApproximateGroundStationPositions, R"doc(
        
    This function returns the positions of DSN ground stations and VLBI stations.

    The function returns ground station positions for all DSN ground stations (as given by :func:`~tudatpy.dynamics.environment_setup.ground_station.get_approximate_dsn_ground_station_positions`) and VLBI stations (as given by :func:`~tudatpy.dynamics.environment_setup.ground_station.get_vlbi_station_positions`).

    Returns
    -------
    dict[str, numpy.ndarray([3,1])]
        Dictionary mapping station name to positions.
        
        )doc" );

    m.def( "approximate_ground_stations_position", &tss::getCombinedApproximateGroundStationPositions, R"doc(No documentation found.)doc" );

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
}

}  // namespace ground_station
}  // namespace environment_setup
}  // namespace dynamics
}  // namespace tudatpy
