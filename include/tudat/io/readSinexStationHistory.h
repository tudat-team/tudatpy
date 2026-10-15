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
struct SinexSourceRecord {
    std::string source;
    std::size_t line = 0;
    std::string text;
};
struct SinexCivilUtc {
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

//! Zero tokens carry their field's meaning; they are never converted to a date/infinity.
enum class SinexEpochMeaning { specified, unspecified, headerStartDefault, headerEndDefault, openStart, openEnd };
struct SinexEpoch {
    std::string original;
    SinexEpochMeaning meaning = SinexEpochMeaning::unspecified;
    std::optional< SinexCivilUtc > civilUtc;
};
struct SinexSiteKey {
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
struct SinexSolutionKey {
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
struct SinexSiteRecord {
    SinexSiteKey key;
    std::string domesId;
    std::string observationTechnique;
    std::string description;
    std::optional< double > approximateLongitudeRadians;
    std::optional< double > approximateLatitudeRadians;
    std::optional< double > approximateHeightMetres;
    SinexSourceRecord provenance;
};
enum class SinexNormalization { si, unresolvedExtension, unsupportedUnit };
struct SinexParameterRecord {
    std::size_t index = 0;
    std::string type;
    SinexSolutionKey key;
    SinexEpoch referenceEpoch;
    std::string originalUnit;
    std::string constraintCode;
    double originalValue = 0.0;
    double originalStandardDeviation = 0.0;
    SinexNormalization normalization = SinexNormalization::unsupportedUnit;
    std::optional< double > valueSi;
    std::optional< double > standardDeviationSi;
    std::string siUnit;
    SinexSourceRecord provenance;
};
struct SinexObservationSpan {
    SinexSolutionKey key;
    std::string observationTechnique;
    SinexEpoch start;
    SinexEpoch end;
    SinexEpoch mean;
    SinexSourceRecord provenance;
};

//! ITRF extension: position and velocity identifiers remain independent.
struct SinexDiscontinuityRecord {
    SinexSolutionKey key;
    std::string observationTechnique;
    std::string kind;
    SinexEpoch start;
    SinexEpoch end;
    std::string description;
    SinexSourceRecord provenance;
};
enum class SinexEccentricityBasis { une, xyz };
struct SinexEccentricityRecord {
    SinexSolutionKey key;
    std::string observationTechnique;
    SinexEpoch start;
    SinexEpoch end;
    SinexEccentricityBasis basis = SinexEccentricityBasis::une;

    //! Metres in the declared component order; UNE is not silently reordered to ENU.
    std::array< double, 3 > componentsMetres{ { 0.0, 0.0, 0.0 } };
    SinexSourceRecord provenance;
};
enum class SinexEquipmentKind { receiver, antenna };
struct SinexEquipmentRecord {
    SinexEquipmentKind kind = SinexEquipmentKind::receiver;
    SinexSolutionKey key;
    std::string observationTechnique;
    SinexEpoch start;
    SinexEpoch end;

    //! Antenna model includes the source's radome suffix; no model/calibration is inferred.
    std::string equipmentType;
    std::string serialNumber;
    std::string firmwareVersion;
    SinexSourceRecord provenance;
};
struct SinexOpaqueBlock {
    std::string name;
    std::vector< SinexSourceRecord > records;
};
struct SinexStationHistory {
    std::string source;
    std::string originalHeader;
    std::string version;
    std::string agency;
    SinexEpoch creationEpoch;
    SinexEpoch headerStart;
    SinexEpoch headerEnd;
    std::string observationTechnique;
    std::size_t declaredParameterCount = 0;

    //! Neutral supplied metadata; absent information stays empty, never inferred from filenames.
    std::string frame;
    std::string release;
    std::string origin;
    std::string tideSystem;
    std::vector< SinexSourceRecord > comments;

    //! Repeated raw site/point keys can identify different techniques/monuments in combined products.
    std::multimap< SinexSiteKey, SinexSiteRecord > sites;
    std::map< std::size_t, SinexParameterRecord > parameters;
    std::vector< SinexObservationSpan > observationSpans;
    std::vector< SinexDiscontinuityRecord > discontinuities;
    std::vector< SinexEccentricityRecord > eccentricities;
    std::vector< SinexEquipmentRecord > equipmentHistories;

    //! Unknown blocks and covariance declarations are retained without evaluation.
    std::vector< SinexOpaqueBlock > opaqueBlocks;
};
struct SinexCartesianSolution {
    SinexSolutionKey key;
    SinexSiteRecord site;
    SinexEpoch positionReferenceEpoch;
    std::array< double, 3 > positionMetres{ { 0.0, 0.0, 0.0 } };
    std::array< std::size_t, 3 > positionParameterIndices{ { 0, 0, 0 } };
    std::array< double, 3 > positionStandardDeviationsMetres{ { 0.0, 0.0, 0.0 } };
    std::array< std::optional< double >, 3 > velocityMetresPerSecond;
    std::array< std::optional< std::size_t >, 3 > velocityParameterIndices;
    std::array< std::optional< double >, 3 > velocityStandardDeviationsMetresPerSecond;
    std::array< std::optional< SinexEpoch >, 3 > velocityReferenceEpochs;
};

//! Strict SINEX 2.02 fixed-width station-product input, including parameter-only IGN extensions.
SinexStationHistory readSinexStationHistory( std::istream& stream, const std::string& source = "<stream>" );
SinexStationHistory readSinexStationHistory( const std::string& fileName );

//! Validate mutable preloaded records without constructing any environment model.
void validateSinexStationHistory( const SinexStationHistory& history );

//! Explicit solution selection; complete positions are required only by this consuming query.
SinexCartesianSolution getSinexCartesianSolution( const SinexStationHistory& history,
                                                  const SinexSolutionKey& key,
                                                  bool requireCompleteVelocity = false );
std::vector< SinexSolutionKey > getSinexCartesianSolutionKeys( const SinexStationHistory& history );

}  // namespace input_output
}  // namespace tudat
#endif
