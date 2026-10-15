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

#include <istream>
#include <map>
#include <utility>
#include <string>
#include <vector>

#include "tudat/basics/timeType.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace input_output
{

//! Hardware identifiers and wavelength of one C0 system configuration.
struct CrdSystemConfigurationRecord {

    //! Identifier referenced by range, calibration and supplement records.
    std::string systemConfigurationId_;

    //! Transmit wavelength in nanometres and component configuration identifiers.
    double transmitWavelengthNm_ = TUDAT_NAN;

    std::vector< std::string > componentConfigurationIds_;
};

//! Station, target and session metadata; applied flags retain -1 for unavailable values.
struct CrdPassConfigurationData {

    //! Source and CRD version, station name, pad, system, occupancy and UTC reference code.
    std::string source_;

    int formatVersion_ = -1;

    std::string stationName_;

    int cdpPadId_ = -1;

    int systemNumber_ = -1;

    int occupancy_ = -1;

    int stationEpochTimeScale_ = -1;

    std::string stationNetwork_;

    //! Target identifiers and data type (0 full rate, 1 normal points, 2 engineering).
    std::string targetName_;

    std::string ilrsSatelliteId_;

    int dataType_ = -1;

    //! Session start and end calendar fields in the declared station UTC reference.
    int startYear_ = -1;

    int startMonth_ = -1;

    int startDay_ = -1;

    int startHour_ = -1;

    int startMinute_ = -1;

    long double startSecond_ = TUDAT_NAN;

    int endYear_ = -1;

    int endMonth_ = -1;

    int endDay_ = -1;

    int endHour_ = -1;

    int endMinute_ = -1;

    long double endSecond_ = TUDAT_NAN;

    //! Release and correction flags from H4; 0 means absent and 1 means already applied.
    int dataRelease_ = -1;

    int troposphericCorrectionApplied_ = -1;

    int centerOfMassCorrectionApplied_ = -1;

    int receiveAmplitudeCorrectionApplied_ = -1;

    int stationSystemDelayApplied_ = -1;

    int spacecraftSystemDelayApplied_ = -1;

    //! H4 range type (2 is two-way) and quality alert (0 is nominal).
    int rangeType_ = -1;

    int dataQualityAlert_ = -1;

    //! C0 configurations, indexed by their published identifiers rather than by station.
    std::map< std::string, CrdSystemConfigurationRecord > systemConfigurations_;
};

//! A record11 normal point; time of flight is in seconds and bin statistics are in picoseconds.
struct CrdNormalPointRecord {

    //! Original seconds of day, measured two-way interval in seconds and half-range in metres.
    long double secondOfDay_ = TUDAT_NAN;

    double twoWayTimeOfFlight_ = TUDAT_NAN;

    double oneWayRange_ = TUDAT_NAN;

    //! The system identifier and epoch event associate each range with its header/configuration.
    std::string systemConfigurationId_ = "";

    int epochEvent_ = -1;

    //! Window length in seconds, number of returns and raw two-way bin RMS in ps.
    double normalPointWindowLength_ = TUDAT_NAN;

    int numberOfReturns_ = -1;

    double binRms_ = TUDAT_NAN;

    //! Shape statistics, peak-minus-mean in ps, return rate in percent and optional signal-to-noise ratio.
    double binSkew_ = TUDAT_NAN;

    double binKurtosis_ = TUDAT_NAN;

    double binPeakMinusMean_ = TUDAT_NAN;

    double returnRate_ = TUDAT_NAN;

    int detectorChannel_ = -1;

    double signalToNoiseRatio_ = TUDAT_NAN;
};

//! A record10 full-rate measurement, retaining the published epoch-event and detector identifiers.
struct CrdFullRateRecord {

    //! Original seconds of day, measured interval in seconds and half-range in metres.
    long double secondOfDay_ = TUDAT_NAN;

    double twoWayTimeOfFlight_ = TUDAT_NAN;

    double oneWayRange_ = TUDAT_NAN;

    //! The system identifier and epoch event associate each range with its header/configuration.
    std::string systemConfigurationId_ = "";

    int epochEvent_ = -1;

    //! Filter, detector and stop identifiers describe the original full-rate return.
    int filterFlag_ = -1;

    int detectorChannel_ = -1;

    int stopNumber_ = -1;
};

//! A record20 weather sample: pressure in hPa, temperature in kelvin and relative humidity in percent.
struct CrdMeteoRecord {

    //! Sample epoch and the recorded atmospheric quantities.
    long double secondOfDay_ = TUDAT_NAN;

    double pressure_ = TUDAT_NAN;

    double temperature_ = TUDAT_NAN;

    double humidity_ = TUDAT_NAN;

    //! 0 denotes a measured weather sample; 1 denotes an interpolated sample.
    int originOfValues_ = -1;
};

//! Supplementary corrections from record12; the H4 flags determine whether corrections were applied.
struct CrdRangeSupplementRecord {

    //! Seconds of day and associated C0 system identifier.
    long double secondOfDay_ = TUDAT_NAN;

    std::string systemConfigurationId_;

    //! One-way atmospheric correction in ps and centre-of-mass correction in metres.
    double troposphericCorrectionPs_ = TUDAT_NAN;

    double centerOfMassCorrection_ = TUDAT_NAN;

    //! Filter value, applied epoch bias in seconds and optional TOF rate in seconds/second.
    double neutralDensityFilter_ = TUDAT_NAN;

    double timeBiasApplied_ = TUDAT_NAN;

    double rangeRate_ = TUDAT_NAN;
};

//! Calibration summary40 or detail41; delays and distribution statistics remain in printed ps units.
struct CrdCalibrationRecord {

    //! Original record type, seconds of day, ranging type and C0 configuration identifier.
    int recordType_ = -1;

    long double secondOfDay_ = TUDAT_NAN;

    int dataType_ = -1;

    std::string systemConfigurationId_;

    //! Recorded/used counts and one-way calibration target distance in metres.
    int numberRecorded_ = -1;

    int numberUsed_ = -1;

    double targetDistance_ = TUDAT_NAN;

    //! System delay, delay shift, RMS, skew, kurtosis and peak-minus-mean.
    double systemDelayPs_ = TUDAT_NAN;

    double delayShiftPs_ = TUDAT_NAN;

    double rmsPs_ = TUDAT_NAN;

    double skew_ = TUDAT_NAN;

    double kurtosis_ = TUDAT_NAN;

    double peakMinusMeanPs_ = TUDAT_NAN;

    //! Calibration type/shift type, detector channel, span and return rate in percent.
    int calibrationType_ = -1;

    int shiftType_ = -1;

    int detectorChannel_ = -1;

    int span_ = -1;

    double returnRate_ = TUDAT_NAN;
};

//! Ordered session records; repeated observation epochs are retained.
struct CrdPassData {
    std::vector< CrdFullRateRecord > fullRateData_;

    std::vector< CrdNormalPointRecord > normalPointData_;

    std::vector< CrdMeteoRecord > meteorologicalData_;

    std::vector< CrdRangeSupplementRecord > rangeSupplementData_;

    std::vector< CrdCalibrationRecord > calibrationData_;
};

//! One CRD session with its station, target, system and correction metadata.
struct CrdPass {
    CrdPassConfigurationData configuration_;

    CrdPassData data_;
};

//! Convert a two-way light time in seconds to one-way range in metres.
double convertCrdTwoWayTimeOfFlightToSlrRange( const double twoWayTimeOfFlight );

//! Read CRDv1/v2 records without applying time, atmospheric, calibration or target corrections.
std::vector< CrdPass > readCrdFile( std::istream& stream, const std::string& source = "<stream>" );

//! Read a CRD file, retaining the filename in each session's source metadata.
std::vector< CrdPass > readCrdFile( const std::string& fileName );

//! Concatenate sessions from the supplied CRD files in input order.
std::vector< CrdPass > readCrdFiles( const std::vector< std::string >& fileNames );

//! Group sessions by their published target names.
std::map< std::string, std::vector< CrdPass > > groupCrdDataPerTarget( const std::vector< CrdPass >& crdPasses );

//! Group sessions using an explicit station-pad to user station-name map.
std::map< std::string, std::vector< CrdPass > > groupCrdDataPerStation(
        const std::vector< CrdPass >& crdPasses,
        const std::map< int, std::string >& monumentIdToGroundStationNameMap );

//! Extract ordered (UTC-coordinate Time since J2000, one-way range in metres) pairs, preserving duplicates.
//! Only declared two-way data and station UTC reference codes 3/4/7 are supported; epoch events are unchanged.
std::vector< std::pair< Time, double > > extractNormalPointMeasurements( const CrdPass& passData );

//! Extract ordered normal points from multiple sessions, without merging equal epochs.
std::vector< std::pair< Time, double > > extractNormalPointMeasurements( const std::vector< CrdPass >& passData );

//! Extract ordered two-way full-rate half-ranges at the original epoch events in the declared station UTC reference.
std::vector< std::pair< Time, double > > extractFullRateMeasurements( const CrdPass& passData );

//! Extract full-rate records from multiple sessions, without merging equal epochs.
std::vector< std::pair< Time, double > > extractFullRateMeasurements( const std::vector< CrdPass >& passData );

//! Return station wavelengths in metres, rejecting stations whose C0 configurations use different wavelengths.
std::map< std::string, double > getStationWavelengths( const std::map< std::string, std::vector< CrdPass > >& groupedData );

}  // namespace input_output

}  // namespace tudat

#endif  // TUDAT_READCRDFILE_H
