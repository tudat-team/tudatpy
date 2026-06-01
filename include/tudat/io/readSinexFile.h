/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_READSINEXFILE_H
#define TUDAT_READSINEXFILE_H

#include <map>
#include <string>

#include <Eigen/Core>

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/basic_astro/timeConversions.h"

namespace tudat
{

namespace input_output
{

struct SinexStationState {
    std::string siteCode_ = "";
    std::string domesId_ = "";
    Eigen::Vector3d position_ = Eigen::Vector3d::Constant( TUDAT_NAN );
    Eigen::Vector3d velocity_ = Eigen::Vector3d::Constant( TUDAT_NAN );
    double referenceEpoch_ = TUDAT_NAN;
};

struct SinexStationEccentricity {
    std::string domesId_ = "";
    int stationCode_ = -1;
    Eigen::Vector3d eccentricity_ = Eigen::Vector3d::Zero( );
    double startEpoch_ = TUDAT_NAN;
    double endEpoch_ = TUDAT_NAN;
    bool hasOpenEnd_ = false;
};

struct IlrsStationRegistryEntry {
    int stationCode_ = -1;
    std::string stationName_ = "";
    std::string domesId_ = "";
    double approximateLongitude_ = TUDAT_NAN;
    double approximateLatitude_ = TUDAT_NAN;
    double approximateHeight_ = TUDAT_NAN;
};

double convertSinexDateTimeToSecondsSinceEpoch( const std::string& dateTime,
                                                const double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

std::map< std::string, SinexStationState > readSinexStationData(
        const std::string& fileName,
        const double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

std::map< std::string, std::vector< SinexStationEccentricity > > readSinexStationEccentricities(
        const std::string& fileName,
        const double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

std::map< int, IlrsStationRegistryEntry > readIlrsStationRegistryFromSinexSiteId( const std::string& fileName );

std::map< std::string, std::string > readDomesIdNumbers( const std::string& fileName );

std::map< int, std::string > readMonumentNumbers( const std::string& fileName );

std::map< std::string, std::string > readGroundStationNames( const std::string& fileName );

}  // namespace input_output

}  // namespace tudat

#endif  // TUDAT_READSINEXFILE_H
