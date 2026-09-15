/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DEFAULTGROUNDSTATIONSETTINGS_H
#define TUDAT_DEFAULTGROUNDSTATIONSETTINGS_H

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

namespace tudat
{

namespace simulation_setup
{

class GroundStationSettings;

/*!
 * Returns a map with the approximate positions of the DSN ground stations, having as key the ground station names. The
 * ground stations are named "DSS-id". The ground station positions are selected according to table 2 of DSN 810-005,
 * 301 Coverage and Geometry, Revision O (2024), DSN/JPL. The positions of the ground stations are specified at 2003.0
 * with respect to ITRF93.
 *
 * @return Map with the ground station positions
 */
const std::map< std::string, Eigen::Vector3d >& getApproximateDsnGroundStationPositions( );

std::map< std::string, Eigen::Vector3d > getCombinedApproximateGroundStationPositions( );

/*!
* Returns the complex id for a DSN ground station from the station name.
* @param stationName Name of the DSN ground station
* @return Complex id for the DSN ground station, this is 10 for Goldstone, 40 for Canberra and 60 for Madrid.

*/
int getDsnComplexId( const std::string& stationName );

/*!
 * Returns the approximate position of the specified ground station. Currently only implemented for DSN stations.
 *
 * @param stationName Station name
 * @return Ground station position
 */
Eigen::Vector3d getApproximateGroundStationPosition( const std::string& stationName );

const std::map< std::string, Eigen::Vector3d >& getVlbiStationPositions( );

const std::map< std::string, Eigen::Vector3d >& getVlbiStationVelocities( );

/*!
 * Returns the velocity for a DSN ground station. The velocities are specified according to table 3 of DSN 810-005,
 * 301 Coverage and Geometry, Revision O (2024), DSN/JPL.
 *
 * @return Velocity for respective station.
 */
Eigen::Vector3d getDsnStationVelocity( const std::string& stationName );

/*!
 * Returns the setting for a DSN ground station. The settings are specified according to table 2 and 3 of DSN 810-005,
 * 301 Coverage and Geometry, Revision O (2024), DSN/JPL. The positions of the ground stations are specified with respect
 * to ITRF2014 and account for their linear motion.
 *
 * @return Ground station settings for respective station.
 */
std::shared_ptr< GroundStationSettings > getDsnStationSetting( const std::string& stationName );

/*!
 * Returns the settings for DSN ground stations. The settings are specified according to table 2 and 3 of DSN 810-005,
 * 301 Coverage and Geometry, Revision O (2024), DSN/JPL. The positions of the ground stations are specified with respect
 * to ITRF2014 and account for their linear motion.
 *
 * @return Vector of ground station settings.
 */
std::vector< std::shared_ptr< GroundStationSettings > > getDsnStationSettings( );

std::vector< std::shared_ptr< GroundStationSettings > > getEvnStationSettings( );

std::vector< std::shared_ptr< GroundStationSettings > > getMPCStationSettings( );

std::vector< std::shared_ptr< GroundStationSettings > > getRadioTelescopeStationSettings( );

std::string getDefaultIlrsSinexStateFilePath( );

std::string getDefaultIlrsSinexEccentricityFilePath( );

std::vector< std::shared_ptr< GroundStationSettings > > getIlrsStationSettingsFromSinexDomes(
        const std::vector< std::string >& domesIds,
        const std::string& sinexStateFile,
        const std::string& sinexEccentricityFile,
        const bool throwExceptionOnMissingData = false );

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_DEFAULTGROUNDSTATIONSETTINGS_H
