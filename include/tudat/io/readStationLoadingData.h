#ifndef TUDAT_READ_STATION_LOADING_DATA_H
#define TUDAT_READ_STATION_LOADING_DATA_H

#include <array>
#include <cstddef>
#include <istream>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace tudat
{
namespace input_output
{

enum class StationLoadingTimeScale { tai, tt, utc, tdb };
enum class StationLoadingBasis { unspecified, geocentricRadialEastNorth, cartesian };
enum class StationLoadingComponent { atmosphere, hydrology, ocean };
enum class StationLoadingSignalPresence { unknown, absent, present };
enum class StationLoadingOrigin { unspecified, wholeEarthCentreOfMass, solidEarthCentreOfMass, centreOfFigure, custom };

//! Authoritative source epoch, retained without a single large floating-point Julian date.
struct StationLoadingEpoch {
    int modifiedJulianDay = 0;
    long double secondsOfDay = 0.0L;
};

//! Scientific inventory is caller supplied; source labels and filenames do not establish it.
struct StationLoadingMetadata {
    std::string source;
    std::vector< std::string > sourceHeaders;
    std::string product;
    std::string release;
    std::string sourceProfile;
    std::string referenceFrame;
    StationLoadingOrigin referenceOrigin = StationLoadingOrigin::unspecified;
    std::string originDefinition;
    std::string oceanResponse;
    std::vector< StationLoadingComponent > components;
    StationLoadingSignalPresence diurnal = StationLoadingSignalPresence::unknown;
    StationLoadingSignalPresence semidiurnal = StationLoadingSignalPresence::unknown;
    StationLoadingSignalPresence annual = StationLoadingSignalPresence::unknown;
    StationLoadingSignalPresence semiannual = StationLoadingSignalPresence::unknown;
    StationLoadingSignalPresence mean = StationLoadingSignalPresence::unknown;
    StationLoadingSignalPresence trend = StationLoadingSignalPresence::unknown;
    std::vector< std::string > removedHarmonics;
    std::string originalDisplacementUnits = "m";
};

struct StationLoadingSample {
    //! Metres, in the data object's explicitly declared source basis/order.
    std::array< double, 3 > displacement{};
    std::size_t sourceLine = 0;
    std::string sourceRecord;

    //! EPHEDISP columns 10--43: preserved as text, never used as the sample epoch.
    std::string informationalEpoch;
};

struct StationLoadingSite {
    //! Internal source label, without an implicit identity mapping to a station name.
    std::string identifier;
    std::array< double, 3 > cartesianPosition{};
    double validityRadiusMetres = 0.0;
    std::size_t sourceLine = 0;
    std::string sourceRecord;
    std::string informationalCoordinates;

    //! Zero-based indexes into epochs. Missing leading/trailing samples or all samples are allowed.
    std::map< std::size_t, StationLoadingSample > samples;
};

//! A generic preloaded series supports an explicit irregular grid as well as EPHEDISP's regular grid.
//! Every site's nonempty sample index range must be contiguous; no internal gap is silently filled.
struct StationLoadingData {
    StationLoadingMetadata metadata;
    StationLoadingTimeScale timeScale = StationLoadingTimeScale::tai;
    StationLoadingBasis basis = StationLoadingBasis::unspecified;
    std::vector< StationLoadingEpoch > epochs;
    std::optional< StationLoadingEpoch > declaredBeginEpoch;
    std::optional< StationLoadingEpoch > declaredEndEpoch;
    std::optional< long double > nominalSampleIntervalSeconds;
    std::map< std::string, StationLoadingSite > sites;
    std::string informationalBeginEpoch;
    std::string informationalEndEpoch;
};

void validateStationLoadingData( const StationLoadingData& data );

//! Strict fixed-column EPHEDISP 2005.06.30 reader; authoritative epochs/indexes are TAI.
//! Both specification T/A and observed IMLS A/T section ordering are supported.
StationLoadingData readEphedispStationLoadingData( std::istream& stream, const StationLoadingMetadata& metadata = {} );
StationLoadingData readEphedispStationLoadingData( const std::string& path, const StationLoadingMetadata& metadata = {} );

//! Explicit documented source profile. Release remains a mandatory caller identifier.
//! The source's whole-Earth CM definition does not assert CE/CF equivalence.
StationLoadingMetadata imlsMerra2AtmosphericLoadingMetadata( const std::string& release );

}  // namespace input_output
}  // namespace tudat
#endif
