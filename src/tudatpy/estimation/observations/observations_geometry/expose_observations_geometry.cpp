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
#include "expose_observations_geometry.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;

namespace tudat
{

namespace simulation_setup
{

std::pair< std::vector< double >, std::vector< Eigen::VectorXd > > getTargetAnglesAndRangeVector(
        const simulation_setup::SystemOfBodies &bodies,
        const std::pair< std::string, std::string > groundStationId,
        const std::string &targetBody,
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

 observation_times : list[astro.time_representation.Time]
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

 observation_times : list[astro.time_representation.Time]
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
