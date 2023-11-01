/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DEFAULTBODIES_H
#define TUDAT_DEFAULTBODIES_H

#include "tudat/simulation/environment_setup/createBodies.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create default settings for a body's atmosphere model.
/*!
 *  Function to create default settings for a body's atmosphere model. Currently set to no
 *  atmosphere, except for Earth, for which a tabulated version of the 1976 Standard Atmosphere is
 *  set.
 *  \param bodyName Name of body for which default atmosphere settings are to be retrieved.
 *  \param initialTime Start time at which environment models in body are to be created
 *  (not currently used by this function, but included for consistency).
 *  \param finalTime End time up to which environment models in body are to be created
 *  (not currently used by this function, but included for consistency).
 *  \return Default settings for a body's atmosphere model.
 */
std::shared_ptr< AtmosphereSettings > getDefaultAtmosphereModelSettings(
        const std::string& bodyName,
        const double initialTime,
        const double finalTime );

std::shared_ptr<RadiationSourceModelSettings> getKnockeEarthRadiationPressureSettings( );

//! Function to create default settings for a body's radiation source model.
/*!
 *  Function to create default settings for a body's radiation source model. Currently set to no
 *  source, except for the Sun.
 *  \param bodyName Name of body for which default radiation source model settings are to be retrieved.
 *  \param initialTime Start time at which environment models in body are to be created
 *  (not currently used by this function, but included for consistency).
 *  \param finalTime End time up to which environment models in body are to be created
 *  (not currently used by this function, but included for consistency).
 *  \return Default settings for a body's radiation source model.
 */
std::shared_ptr< RadiationSourceModelSettings > getDefaultRadiationSourceModelSettings(
        const std::string& bodyName,
        const double initialTime,
        const double finalTime );

//! Function to create default settings for a body's ephemeris.
/*!
 *  Function to create default settings for a body's ephemeris without a limitation on the time interval.
 *  \param bodyName Name of body for which default ephemeris settings are to be retrieved.
 *  \return Default settings for a body's ephemeris.
 */
std::shared_ptr< EphemerisSettings > getDefaultEphemerisSettings(
        const std::string& bodyName,
        const std::string& baseFrameOrientation = "ECLIPJ2000",
        const std::string& originatingNameBodyName = "" );

//! Function to create default settings for a body's ephemeris.
/*!
 *  Function to create default settings for a body's ephemeris. Currently set to a
 *  creating a 6th order Lagrange interpolator from Spice.
 *  \param bodyName Name of body for which default ephemeris settings are to be retrieved.
 *  \param initialTime Start time at which ephemeris is to be created.
 *  \param finalTime End time up to which ephemeris is to be created.
 *  \param timeStep Time step with which interpolated data from Spice should be created.
 *  \return Default settings for a body's ephemeris.
 */
std::shared_ptr< EphemerisSettings > getDefaultEphemerisSettings(
        const std::string& bodyName,
        const double initialTime,
        const double finalTime,
        const std::string& baseFrameOrientation = "ECLIPJ2000",
        const std::string& originatingNameBodyName = "",
        const double timeStep = 300.0 );

//! Function to create default settings for a body's gravity field model.
/*!
 *  Function to create default settings for a body's gravity field model. Currently set to
 *  a point mass gravty field, with the gravitational parameter obtained from Spice.
 *  \param bodyName Name of body for which default gravity field settings are to be retrieved.
 *  \param initialTime Start time at which environment models in body are to be created
 *  (not currently used by this function, but included for consistency).
 *  \param finalTime End time up to which environment models in body are to be created
 *  (not currently used by this function, but included for consistency).
 *  \return Default settings for a body's gravity field model.
 */
std::shared_ptr< GravityFieldSettings > getDefaultGravityFieldSettings(
        const std::string& bodyName,
        const double initialTime,
        const double finalTime );

//! Function to create default settings for a body's rotation model.
/*!
 *  Function to create default settings for a body's rotation model. Currently set to
 *  a rotation model taken directly from Spice
 *  \param bodyName Name of body for which default rotation model settings are to be retrieved.
 *  \param initialTime Start time at which environment models in body are to be created
 *  (not currently used by this function, but included for consistency).
 *  \param finalTime End time up to which environment models in body are to be created
 *  (not currently used by this function, but included for consistency).
 *  \return Default settings for a body's rotation model.
 */
std::shared_ptr< RotationModelSettings > getDefaultRotationModelSettings(
        const std::string& bodyName,
        const double initialTime,
        const double finalTime,
        const std::string& baseFrameOrientation = "ECLIPJ2000",
        const std::string& spiceBodyName = "" );

double marsTimeDependentPhaseAngleCorrectionFunction( const double secondsSinceJ2000 );

std::shared_ptr< RotationModelSettings > getHighAccuracyMarsRotationModel(
    const std::string& baseFrameOrientation = "ECLIPJ2000",
    const std::string& targetFrameOrientation = "Mars_Fixed" );

//! Function to create default settings for a body's shape model.
/*!
 *  Function to create default settings for a body's shape model. Currently set to
 *  a spherical model, with the radius taken from Spice
 *  \param bodyName Name of body for which default shape model settings are to be retrieved.
 *  \param initialTime Start time at which environment models in body are to be created
 *  (not currently used by this function, but included for consistency).
 *  \param finalTime End time up to which environment models in body are to be created
 *  (not currently used by this function, but included for consistency).
 */
std::shared_ptr< BodyShapeSettings > getDefaultBodyShapeSettings(
        const std::string& bodyName,
        const double initialTime,
        const double finalTime );

//! Function to create default settings for a single for body.
/*!
 *  Function to create default settings for a single body from which to create a object using
 *  the code in createBodies.h/.cpp. This function is included to streamline and simplify the
 *  creation of typical celestial bodies. The default settings for the various
 *  environment models of the body are defined in the various functions defined in this file.
 *  \param body Name of body for which default settings are to be retrieved.
 *  \param initialTime Start time at which environment models in body are to be created
 *  (included as some environment models require e.g., interpolators to be created over
 *  a certain time period).
 *  \param finalTime End time up to which environment models in body are to be created
 *  (included as some environment models require e.g., interpolators to be created over
 *  a certain time period).
 *  \param timeStep Time step with which interpolated data from Spice should be created.
 */
std::shared_ptr< BodySettings > getDefaultSingleBodySettings(
        const std::string& body,
        const double initialTime,
        const double finalTime,
        const std::string& baseFrameOrientation = "ECLIPJ2000",
        const double timeStep = 300.0 );

std::shared_ptr< BodySettings > getDefaultSingleAlternateNameBodySettings(
    const std::string& body,
    const std::string& originatingName,
    const double initialTime,
    const double finalTime,
    const std::string& baseFrameOrientation = "ECLIPJ2000",
    const double timeStep = 300.0 );

std::shared_ptr< BodySettings > getDefaultSingleBodySettings(
        const std::string& bodyName,
        const std::string& baseFrameOrientation = "ECLIPJ2000" );

std::shared_ptr< BodySettings > getDefaultSingleAlternateNameBodySettings(
    const std::string& bodyName,
    const std::string& originatingName,
    const std::string& baseFrameOrientation = "ECLIPJ2000" );

//! Function to create default settings from which to create a set of body objects.
/*!
 *  Function to create default settings from which to create a set of body objects using
 *  the code in createBodies.h/.cpp. This function is included to streamline and simplify the
 *  creation of typical celestial bodies. The default settings for the various
 *  environment models of the body are defined in the various functions defined in this file.
 *  \param bodies List of bodies for which default settings are to be retrieved.
 *  \param initialTime Start time at which environment models in body are to be created
 *  (included as some environment models require e.g., interpolators to be created over
 *  a certain time period).
 *  \param finalTime End time up to which environment models in body are to be created
 *  (included as some environment models require e.g., interpolators to be created over
 *  a certain time period).
 *  \param timeStep Time step with which interpolated data from Spice should be created.
 *  \return Default settings from which to create a set of body objects.
 */
BodyListSettings getDefaultBodySettings(
        const std::vector< std::string >& bodies,
        const double initialTime,
        const double finalTime,
        const std::string baseFrameOrigin = "SSB",
        const std::string baseFrameOrientation = "ECLIPJ2000",
        const double timeStep = 300.0 );

//! Function to create default settings from which to create a set of body objects, without stringent limitations on
//! time-interval of validity of environment.
/*!
 *  Function to create default settings from which to create a set of body objects using
 *  the code in createBodies.h/.cpp. This function is included to streamline and simplify the
 *  creation of typical celestial bodies. The default settings for the various
 *  environment models of the body are defined in the various functions defined in this file.
 *  \param bodies List of bodies for which default settings are to be retrieved.
 *  \return Default settings from which to create a set of body objects.
 */
BodyListSettings getDefaultBodySettings(
        const std::vector< std::string >& bodies,
        const std::string baseFrameOrigin = "SSB",
        const std::string baseFrameOrientation = "ECLIPJ2000" );

/*!
 * Returns a map with the approximate positions of the DSN ground stations, having as key the ground station names. The
 * ground stations are named "DSS-id". The ground station positions are selected according to table 2 of DSN 810-005,
 * 301 Coverage and Geometry, Revision K (2016), DSN/JPL. The positions of the ground stations are specified at 2003.0
 * with respect to ITRF93.
 *
 * @return Map with the ground station positions
 */
std::map< std::string, Eigen::Vector3d > getApproximateDsnGroundStationPositions( );

/*!
 * Returns the default DSN station names per DSN station complex id. Stations are named as "DSS-i", following the
 * nomenclature used when retrieving the default DSN ground station settings.
 */
inline std::map< int, std::vector< std::string > > getDefaultDsnStationNamesPerComplex( )
{
    std::map< int, std::vector< std::string > > stationsPerComplex;
    stationsPerComplex[ 10 ] = { "DSS-13", "DSS-14", "DSS-15", "DSS-24", "DSS-25", "DSS-26", "DSS-27" };
    stationsPerComplex[ 40 ] = { "DSS-34", "DSS-35", "DSS-36", "DSS-43", "DSS-45" };
    stationsPerComplex[ 60 ] = { "DSS-54", "DSS-55", "DSS-63", "DSS-65" };

    return stationsPerComplex;
}

/*!
 * Returns the approximate position of the specified ground station. Currently only implemented for DSN stations.
 *
 * @param stationName Station name
 * @return Ground station position.
 */
Eigen::Vector3d getApproximateGroundStationPosition( std::string stationName );

/*!
 * Returns the settings for DSN ground stations. The settings are specified according to table 2 and 3 of DSN 810-005,
 * 301 Coverage and Geometry, Revision K (2016), DSN/JPL. The positions of the ground stations are specified with respect
 * to ITRF2014 and account for their linear motion.
 *
 * @return Vector of ground station settings.
 */
std::vector< std::shared_ptr< GroundStationSettings > > getDsnStationSettings( );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_DEFAULTBODIES_H
