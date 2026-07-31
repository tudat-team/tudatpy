/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_READCRDFILE_H
#define TUDAT_READCRDFILE_H

#include <map>
#include <string>
#include <vector>

#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace input_output
{

struct CrdPassConfigurationData {
    std::string stationName_ = "";
    int cdpPadId_ = -1;
    std::string targetName_ = "";
    int startYear_ = -1;
    int startMonth_ = -1;
    int startDay_ = -1;
    int endYear_ = -1;
    int endMonth_ = -1;
    int endDay_ = -1;
    double transmitWavelengthNm_ = TUDAT_NAN;
};

struct CrdNormalPointRecord {
    double secondOfDay_ = TUDAT_NAN;
    double twoWayTimeOfFlight_ = TUDAT_NAN;
    double oneWayRange_ = TUDAT_NAN;
    std::string systemConfigurationId_ = "";
    int epochEvent_ = -1;
    double normalPointWindowLength_ = TUDAT_NAN;
    int numberOfReturns_ = -1;
    double binRms_ = TUDAT_NAN;
};

struct CrdFullRateRecord {
    double secondOfDay_ = TUDAT_NAN;
    double twoWayTimeOfFlight_ = TUDAT_NAN;
    double oneWayRange_ = TUDAT_NAN;
    std::string systemConfigurationId_ = "";
    int epochEvent_ = -1;
    int filterFlag_ = -1;
    int detectorChannel_ = -1;
    int stopNumber_ = -1;
};

struct CrdMeteoRecord {
    double secondOfDay_ = TUDAT_NAN;
    double pressure_ = TUDAT_NAN;
    double temperature_ = TUDAT_NAN;
    double humidity_ = TUDAT_NAN;
};

struct CrdPassData {
    std::vector< CrdFullRateRecord > fullRateData_;
    std::vector< CrdNormalPointRecord > normalPointData_;
    std::vector< CrdMeteoRecord > meteorologicalData_;
};

struct CrdPass {
    CrdPassConfigurationData configuration_;
    CrdPassData data_;
};

double convertCrdTwoWayTimeOfFlightToSlrRange( const double twoWayTimeOfFlight );

std::vector< CrdPass > readCrdFile( const std::string& fileName );

std::vector< CrdPass > readCrdFiles( const std::vector< std::string >& fileNames );

std::map< std::string, std::vector< CrdPass > > groupCrdDataPerTarget( const std::vector< CrdPass >& crdPasses );

std::map< std::string, std::vector< CrdPass > > groupCrdDataPerStation(
        const std::vector< CrdPass >& crdPasses,
        const std::map< int, std::string >& monumentIdToGroundStationNameMap );

std::map< double, double > extractNormalPointMeasurements( const CrdPass& passData );

std::map< double, double > extractNormalPointMeasurements( const std::vector< CrdPass >& passData );

std::map< double, double > extractFullRateMeasurements( const CrdPass& passData );

std::map< double, double > extractFullRateMeasurements( const std::vector< CrdPass >& passData );

std::map< std::string, double > getStationWavelengths( const std::map< std::string, std::vector< CrdPass > >& groupedData );

}  // namespace input_output

}  // namespace tudat

#endif  // TUDAT_READCRDFILE_H
