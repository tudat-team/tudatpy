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
#include "expose_observations_geometry.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/astro/ground_stations/pointingAnglesCalculator.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;

namespace tudat
{

namespace simulation_setup
{

Eigen::Vector2d inertialVectorToRightAscensionDeclination( const Eigen::Vector3d& inertialVector )
{
    Eigen::Vector3d normalizedVector = inertialVector.normalized( );
    double rightAscension = std::atan2( normalizedVector.y( ), normalizedVector.x( ) );
    double declination = std::asin( normalizedVector.z( ) );
    return ( Eigen::Vector2d( ) << rightAscension, declination ).finished( );
}

Eigen::Vector3d rightAscensionDeclinationToInertialUnitVector( const Eigen::Vector2d& rightAscensionDeclination )
{
    const double rightAscension = rightAscensionDeclination.x( );
    const double declination = rightAscensionDeclination.y( );
    return ( Eigen::Vector3d( ) << std::cos( declination ) * std::cos( rightAscension ),
             std::cos( declination ) * std::sin( rightAscension ),
             std::sin( declination ) )
            .finished( );
}

std::shared_ptr< ground_stations::PointingAnglesCalculator > getStationPointingAnglesCalculator(
        const simulation_setup::SystemOfBodies& bodies,
        const std::pair< std::string, std::string >& groundStationId )
{
    return bodies.at( groundStationId.first )->getGroundStation( groundStationId.second )->getPointingAnglesCalculator( );
}

Eigen::Vector2d inertialVectorToAzimuthElevation( const simulation_setup::SystemOfBodies& bodies,
                                                  const std::pair< std::string, std::string > groundStationId,
                                                  const double time,
                                                  const Eigen::Vector3d& inertialVector )
{
    std::shared_ptr< ground_stations::PointingAnglesCalculator > pointingAnglesCalculator =
            getStationPointingAnglesCalculator( bodies, groundStationId );
    std::pair< double, double > elevationAzimuth = pointingAnglesCalculator->calculatePointingAngles( inertialVector, time );
    return ( Eigen::Vector2d( ) << elevationAzimuth.second, elevationAzimuth.first ).finished( );
}

Eigen::Vector3d azimuthElevationToInertialUnitVector( const simulation_setup::SystemOfBodies& bodies,
                                                      const std::pair< std::string, std::string > groundStationId,
                                                      const double time,
                                                      const Eigen::Vector2d& azimuthElevation )
{
    std::shared_ptr< ground_stations::PointingAnglesCalculator > pointingAnglesCalculator =
            getStationPointingAnglesCalculator( bodies, groundStationId );

    const double azimuth = azimuthElevation.x( );
    const double elevation = azimuthElevation.y( );
    Eigen::Vector3d topocentricVector = ( Eigen::Vector3d( ) << std::sin( azimuth ) * std::cos( elevation ),
                                          std::cos( azimuth ) * std::cos( elevation ),
                                          std::sin( elevation ) )
                                                .finished( );

    Eigen::Quaterniond rotationFromInertialToTopocentric =
            pointingAnglesCalculator->getRotationFromBodyFixedToTopoCentricFrame( )( time ) *
            pointingAnglesCalculator->getRotsationFromInertialToBodyFixedFrame( )( time );
    return rotationFromInertialToTopocentric.inverse( ) * topocentricVector;
}

Eigen::Vector2d rightAscensionDeclinationToAzimuthElevation( const simulation_setup::SystemOfBodies& bodies,
                                                             const std::pair< std::string, std::string > groundStationId,
                                                             const double time,
                                                             const Eigen::Vector2d& rightAscensionDeclination )
{
    return inertialVectorToAzimuthElevation(
            bodies, groundStationId, time, rightAscensionDeclinationToInertialUnitVector( rightAscensionDeclination ) );
}

Eigen::Vector2d azimuthElevationToRightAscensionDeclination( const simulation_setup::SystemOfBodies& bodies,
                                                             const std::pair< std::string, std::string > groundStationId,
                                                             const double time,
                                                             const Eigen::Vector2d& azimuthElevation )
{
    return inertialVectorToRightAscensionDeclination(
            azimuthElevationToInertialUnitVector( bodies, groundStationId, time, azimuthElevation ) );
}

std::pair< std::vector< double >, std::vector< Eigen::VectorXd > > getTargetAnglesAndRangeVector(
        const simulation_setup::SystemOfBodies& bodies,
        const std::pair< std::string, std::string > groundStationId,
        const std::string& targetBody,
        const std::vector< double > times,
        const bool transmittingToTarget )
{
    std::map< double, Eigen::VectorXd > targetAnglesAndRange =
            getTargetAnglesAndRange( bodies, groundStationId, targetBody, times, transmittingToTarget );
    return std::make_pair( utilities::createVectorFromMapKeys( targetAnglesAndRange ),
                           utilities::createVectorFromMapValues( targetAnglesAndRange ) );
}

}  // namespace simulation_setup

}  // namespace tudat

namespace tudatpy
{
namespace estimation
{
namespace observations
{
namespace observations_geometry
{

void expose_observations_geometry( py::module& m )
{
    m.def( "inertial_vector_to_azimuth_elevation",
           &tss::inertialVectorToAzimuthElevation,
           py::arg( "bodies" ),
           py::arg( "station_id" ),
           py::arg( "time" ),
           py::arg( "inertial_vector" ),
           R"doc(

 Convert an inertial station-to-target vector to local azimuth and elevation at a ground station.

 The supplied vector is rotated to the local topocentric frame of ``station_id`` at ``time``. For topocentric
 direction components :math:`[e,n,u]^T`, the angular convention is

 .. math::
    A &= \operatorname{atan2}(e,n)\\
    E &= \operatorname{atan2}(u,\sqrt{e^2+n^2})

 Parameters
 ----------
 bodies : :class:`~tudatpy.dynamics.environment.SystemOfBodies`
     Environment from which the ground station and frame rotations are retrieved. The target direction is supplied by ``inertial_vector``.

 station_id : tuple[ str, str ]
     Pair ``(body_name, station_name)`` identifying the ground station.

 time : float
     Epoch at which the station frame orientation is evaluated.

 inertial_vector : numpy.ndarray
     Three-element Cartesian vector from the ground station to the target, expressed in the inertial frame.
     Its magnitude is ignored.

 Returns
 -------
 numpy.ndarray
     Two-element vector ``[azimuth, elevation]`` in radians.

     )doc" );

    m.def( "azimuth_elevation_to_inertial_unit_vector",
           &tss::azimuthElevationToInertialUnitVector,
           py::arg( "bodies" ),
           py::arg( "station_id" ),
           py::arg( "time" ),
           py::arg( "azimuth_elevation" ),
           R"doc(

 Convert local azimuth and elevation at a ground station to an inertial station-to-target unit vector.

 For azimuth :math:`A` and elevation :math:`E`, the topocentric unit vector is

 .. math::
    \hat{\mathbf{u}}_{T} =
    [\sin A\cos E,\ \cos A\cos E,\ \sin E]^T

 It is then rotated from the station topocentric frame to the inertial frame.

 Parameters
 ----------
 bodies : :class:`~tudatpy.dynamics.environment.SystemOfBodies`
     Environment from which the ground station and frame rotations are retrieved. The target direction is supplied by ``azimuth_elevation``.

 station_id : tuple[ str, str ]
     Pair ``(body_name, station_name)`` identifying the ground station.

 time : float
     Epoch at which the station frame orientation is evaluated.

 azimuth_elevation : numpy.ndarray
     Two-element vector ``[azimuth, elevation]`` in radians, defining the local direction from the ground station to the target.

 Returns
 -------
 numpy.ndarray
     Three-element Cartesian unit vector from the ground station to the target, expressed in the inertial frame.

     )doc" );

    m.def( "right_ascension_declination_to_azimuth_elevation",
           &tss::rightAscensionDeclinationToAzimuthElevation,
           py::arg( "bodies" ),
           py::arg( "station_id" ),
           py::arg( "time" ),
           py::arg( "right_ascension_declination" ),
           R"doc(

 Convert an inertial right ascension and declination direction to local azimuth and elevation at a ground station.

 For right ascension :math:`\alpha` and declination :math:`\delta`, the inertial station-to-target
 unit vector is :math:`[\cos\delta\cos\alpha,\ \cos\delta\sin\alpha,\ \sin\delta]^T`. This vector is
 then passed to :func:`~tudatpy.estimation.observations.observations_geometry.inertial_vector_to_azimuth_elevation`.

 Parameters
 ----------
 bodies : :class:`~tudatpy.dynamics.environment.SystemOfBodies`
     Environment from which the ground station and frame rotations are retrieved. The target direction is supplied by ``right_ascension_declination``.

 station_id : tuple[ str, str ]
     Pair ``(body_name, station_name)`` identifying the ground station.

 time : float
     Epoch at which the station frame orientation is evaluated.

 right_ascension_declination : numpy.ndarray
     Two-element vector ``[right_ascension, declination]`` in radians, defining the inertial direction from the ground station to the target.

 Returns
 -------
 numpy.ndarray
     Two-element vector ``[azimuth, elevation]`` in radians.

     )doc" );

    m.def( "azimuth_elevation_to_right_ascension_declination",
           &tss::azimuthElevationToRightAscensionDeclination,
           py::arg( "bodies" ),
           py::arg( "station_id" ),
           py::arg( "time" ),
           py::arg( "azimuth_elevation" ),
           R"doc(

 Convert local azimuth and elevation at a ground station to inertial right ascension and declination.

 The local direction is first converted to an inertial unit vector :math:`\hat{\mathbf{u}}` using
 :func:`~tudatpy.estimation.observations.observations_geometry.azimuth_elevation_to_inertial_unit_vector`.
 The returned angles are :math:`\alpha=\operatorname{atan2}(\hat{u}_{y},\hat{u}_{x})` and
 :math:`\delta=\operatorname{asin}(\hat{u}_{z})`.

 Parameters
 ----------
 bodies : :class:`~tudatpy.dynamics.environment.SystemOfBodies`
     Environment from which the ground station and frame rotations are retrieved. The target direction is supplied by ``azimuth_elevation``.

 station_id : tuple[ str, str ]
     Pair ``(body_name, station_name)`` identifying the ground station.

 time : float
     Epoch at which the station frame orientation is evaluated.

 azimuth_elevation : numpy.ndarray
     Two-element vector ``[azimuth, elevation]`` in radians, defining the local direction from the ground station to the target.

 Returns
 -------
 numpy.ndarray
     Two-element vector ``[right_ascension, declination]`` in radians.

     )doc" );

    m.def( "compute_target_angles_and_range_vectors",
           &tss::getTargetAnglesAndRangeVector,
           py::arg( "bodies" ),
           py::arg( "station_id" ),
           py::arg( "target_body" ),
           py::arg( "observation_times" ),
           py::arg( "is_station_transmitting" ),
           R"doc(

 Function to compute the azimuth angle, elevation angle and range at a ground station.

 Function to compute the azimuth angle, elevation angle and range at a ground station. This functions is provided as a function of
 convenience, to prevent users having to manually define the relevant settings for this often-needed functionality. This function
 takes an observing station and a target body as input, and provides the observed angles and current range (without correction for aberrations, with correction for light time)
 as observed at that station


 Parameters
 ----------
 bodies : SystemOfBodies
     System of bodies that defines the full physical environment

 station_id : tuple[ str, str]
     Identifier for the observing station, as a pair of strings: the body name and the station name.

 target_body : str
     Name of body which is observed by ground station

 observation_times : list[:class:`~tudatpy.astro.time_representation.Time`]
     List of times at which the ground station observations are to be analyzed

 is_station_transmitting : bool
     Boolean defining whether the observation times define times at which the station is transmitting to, or receiving from, the ground station.
     This has an impact on the whether the light-time is computed forward or backward in time from the ground station to the target

 Returns
 -------
 dict[astro.time_representation.Time,numpy.ndarray[numpy.float64[3, 1]]]
     Dictionary with the required output. Key defines the observation time, the value is an array of size three containing entry 0 - elevation angle, entry 1 - azimuth angle, entry 2 - range






     )doc" );

    // observation geometry
    m.def( "compute_target_angles_and_range",
           &tss::getTargetAnglesAndRange,
           py::arg( "bodies" ),
           py::arg( "station_id" ),
           py::arg( "target_body" ),
           py::arg( "observation_times" ),
           py::arg( "is_station_transmitting" ),
           R"doc(

 Function to compute the azimuth angle, elevation angle and range at a ground station.

 Function to compute the azimuth angle, elevation angle and range at a ground station. This functions is provided as a function of
 convenience, to prevent users having to manually define the relevant settings for this often-needed functionality. This function
 takes an observing station and a target body as input, and provides the observed angles and current range (without correction for aberrations, with correction for light time)
 as observed at that station


 Parameters
 ----------
 bodies : SystemOfBodies
     System of bodies that defines the full physical environment

 station_id : tuple[ str, str]
     Identifier for the observing station, as a pair of strings: the body name and the station name.

 target_body : str
     Name of body which is observed by ground station

 observation_times : list[:class:`~tudatpy.astro.time_representation.Time`]
     List of times at which the ground station observations are to be analyzed

 is_station_transmitting : bool
     Boolean defining whether the observation times define times at which the station is transmitting to, or receiving from, the ground station.
     This has an impact on the whether the light-time is computed forward or backward in time from the ground station to the target

 Returns
 -------
 dict[float,numpy.ndarray[numpy.float64[3, 1]]]
     Dictionary with the required output. Key defines the observation time, the value is an array of size three containing entry 0 - elevation angle, entry 1 - azimuth angle, entry 2 - range






     )doc" );
}

}  // namespace observations_geometry
}  // namespace observations
}  // namespace estimation
}  // namespace tudatpy
