/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_READ_SINEX_STATION_HISTORY_H
#define TUDAT_READ_SINEX_STATION_HISTORY_H

#include <array>
#include <cstddef>
#include <istream>
#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

namespace tudat
{
namespace input_output
{

//! SINEX epoch (YY:DDD:SSSSS), as calendar year, day of year and second of day in UTC.
struct SinexCivilUtc {
    //! Year (four digits), day of year and second of day.
    int year = 0;
    int dayOfYear = 0;
    int secondsOfDay = 0;

    bool operator==( const SinexCivilUtc& other ) const
    {
        return std::tie( year, dayOfYear, secondsOfDay ) == std::tie( other.year, other.dayOfYear, other.secondsOfDay );
    }
    bool operator<( const SinexCivilUtc& other ) const
    {
        return std::tie( year, dayOfYear, secondsOfDay ) < std::tie( other.year, other.dayOfYear, other.secondsOfDay );
    }
};

//! Meaning of a SINEX epoch field: a date, or one of the meanings of an all-zero epoch (not given, the start or end of the
//! data in the file, or an open start or end of a validity interval). An all-zero epoch is never converted to a date.
enum class SinexEpochMeaning { specified, unspecified, headerStartDefault, headerEndDefault, openStart, openEnd };

//! Epoch field of a SINEX file. An all-zero epoch is not a date; its meaning depends on the field it is in.
struct SinexEpoch {
    //! Whether the field gives a date, or has one of the special meanings of an all-zero epoch.
    SinexEpochMeaning meaning = SinexEpochMeaning::unspecified;

    //! Date and time of the epoch, if the field gives one.
    std::optional< SinexCivilUtc > civilUtc;
};

//! Identifier of a station in a SINEX file: site code and point code.
struct SinexSiteKey {
    //! Four-character site code and two-character point code.
    std::string siteCode;
    std::string pointCode;

    bool operator==( const SinexSiteKey& other ) const
    {
        return siteCode == other.siteCode && pointCode == other.pointCode;
    }
    bool operator<( const SinexSiteKey& other ) const
    {
        return std::tie( siteCode, pointCode ) < std::tie( other.siteCode, other.pointCode );
    }
};

//! Identifier of one solution (one time interval of estimated coordinates) of a station in a SINEX file.
struct SinexSolutionKey {
    //! Station, and solution number of this station.
    SinexSiteKey site;
    std::string solutionId;

    bool operator==( const SinexSolutionKey& other ) const
    {
        return site == other.site && solutionId == other.solutionId;
    }
    bool operator<( const SinexSolutionKey& other ) const
    {
        return std::tie( site.siteCode, site.pointCode, solutionId ) <
                std::tie( other.site.siteCode, other.site.pointCode, other.solutionId );
    }
};

//! Station description from the SITE/ID block of a SINEX file.
struct SinexSiteRecord {
    //! Station identifier, DOMES number, observation technique (e.g. P for GNSS) and free-text description.
    SinexSiteKey key;
    std::string domesId;
    std::string observationTechnique;
    std::string description;

    //! Approximate longitude and latitude (radians) and height (metres) given in the file, if any.
    std::optional< double > approximateLongitudeRadians;
    std::optional< double > approximateLatitudeRadians;
    std::optional< double > approximateHeightMetres;
};

//! Estimated parameter from the SOLUTION/ESTIMATE block of a SINEX file (for instance one coordinate of a station).
struct SinexParameterRecord {
    //! Parameter number in the file, parameter type (e.g. STAX, VELX), station solution and reference epoch.
    std::size_t index = 0;
    std::string type;
    SinexSolutionKey key;
    SinexEpoch referenceEpoch;

    //! Unit, constraint code, value and standard deviation, as given in the file.
    std::string originalUnit;
    std::string constraintCode;
    double originalValue = 0.0;
    double originalStandardDeviation = 0.0;

    //! Value, standard deviation and SI unit, if the value could be converted to SI units.
    std::optional< double > valueSi;
    std::optional< double > standardDeviationSi;
    std::string siUnit;
};

//! Time span of the observations used for one station solution (SOLUTION/EPOCHS block of a SINEX file).
struct SinexObservationSpan {
    //! Station solution and observation technique.
    SinexSolutionKey key;
    std::string observationTechnique;

    //! First, last and mean epoch of the observations.
    SinexEpoch start;
    SinexEpoch end;
    SinexEpoch mean;
};

//! Discontinuity in the coordinates of a station, from the SOLUTION/DISCONTINUITY block (an ITRF extension of SINEX).
//! Position and velocity discontinuities are separate records.
struct SinexDiscontinuityRecord {
    //! Station solution, observation technique and kind of discontinuity (P for position, V for velocity).
    SinexSolutionKey key;
    std::string observationTechnique;
    std::string kind;

    //! Validity interval of the solution, and the reason for the discontinuity.
    SinexEpoch start;
    SinexEpoch end;
    std::string description;
};

//! Component order of a SINEX eccentricity: up, north, east (UNE) or Cartesian X, Y, Z.
enum class SinexEccentricityBasis { une, xyz };

//! Offset of a station's reference point from its marker, from the SITE/ECCENTRICITY block of a SINEX file.
struct SinexEccentricityRecord {
    //! Station solution, observation technique and validity interval of the offset.
    SinexSolutionKey key;
    std::string observationTechnique;
    SinexEpoch start;
    SinexEpoch end;

    //! Component order of the offset.
    SinexEccentricityBasis basis = SinexEccentricityBasis::une;

    //! Metres in the declared component order; UNE is not silently reordered to ENU.
    std::array< double, 3 > componentsMetres{ { 0.0, 0.0, 0.0 } };
};

//! Kind of equipment in a SINEX equipment record.
enum class SinexEquipmentKind { receiver, antenna };

//! Receiver or antenna installed at a station, from the SITE/RECEIVER or SITE/ANTENNA block of a SINEX file.
struct SinexEquipmentRecord {
    //! Kind of equipment, station solution, observation technique and installation interval.
    SinexEquipmentKind kind = SinexEquipmentKind::receiver;
    SinexSolutionKey key;
    std::string observationTechnique;
    SinexEpoch start;
    SinexEpoch end;

    //! Antenna model includes the source's radome suffix; no model/calibration is inferred.
    std::string equipmentType;
    std::string serialNumber;
    std::string firmwareVersion;
};

//! Contents of a SINEX file with station coordinates and station history.
struct SinexStationHistory {
    //! Name of the file (or other source), SINEX version and the agency that made the file.
    std::string source;
    std::string version;
    std::string agency;

    //! Creation epoch of the file, and start and end of the data it contains.
    SinexEpoch creationEpoch;
    SinexEpoch headerStart;
    SinexEpoch headerEnd;

    //! Observation technique(s) and number of estimated parameters, as given in the header.
    std::string observationTechnique;
    std::size_t declaredParameterCount = 0;

    //! Neutral supplied metadata; absent information stays empty, never inferred from filenames.
    std::string frame;
    std::string release;
    std::string origin;
    std::string tideSystem;

    //! Repeated raw site/point keys can identify different techniques/monuments in combined products.
    std::multimap< SinexSiteKey, SinexSiteRecord > sites;
    std::map< std::size_t, SinexParameterRecord > parameters;
    std::vector< SinexObservationSpan > observationSpans;
    std::vector< SinexDiscontinuityRecord > discontinuities;
    std::vector< SinexEccentricityRecord > eccentricities;
    std::vector< SinexEquipmentRecord > equipmentHistories;
};

//! Cartesian position (and, if available, velocity) of one station solution, collected from a SINEX file.
struct SinexCartesianSolution {
    //! Station solution and its description.
    SinexSolutionKey key;
    SinexSiteRecord site;

    //! Reference epoch, position (metres), parameter numbers in the file and standard deviations (metres) of the position.
    SinexEpoch positionReferenceEpoch;
    std::array< double, 3 > positionMetres{ { 0.0, 0.0, 0.0 } };
    std::array< std::size_t, 3 > positionParameterIndices{ { 0, 0, 0 } };
    std::array< double, 3 > positionStandardDeviationsMetres{ { 0.0, 0.0, 0.0 } };

    //! Velocity (metres per second), parameter numbers, standard deviations and reference epochs per component, if given.
    std::array< std::optional< double >, 3 > velocityMetresPerSecond;
    std::array< std::optional< std::size_t >, 3 > velocityParameterIndices;
    std::array< std::optional< double >, 3 > velocityStandardDeviationsMetresPerSecond;
    std::array< std::optional< SinexEpoch >, 3 > velocityReferenceEpochs;
};

//! Read the station blocks of a SINEX (2.0x) file, including the IGN parameter extensions; other blocks are ignored.
SinexStationHistory readSinexStationHistory( std::istream& stream, const std::string& source = "<stream>" );
SinexStationHistory readSinexStationHistory( const std::string& fileName );

//! Collect the Cartesian position (and velocity) of one station solution; the position must be complete.
SinexCartesianSolution getSinexStationCoordinateSolution( const SinexStationHistory& history,
                                                          const SinexSolutionKey& key,
                                                          bool requireCompleteVelocity = false );
std::vector< SinexSolutionKey > getSinexStationCoordinateSolutionKeys( const SinexStationHistory& history );

}  // namespace input_output
}  // namespace tudat
#endif
