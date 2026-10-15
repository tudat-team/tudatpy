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

//! Time scale of the epochs of a station loading series.
enum class StationLoadingTimeScale { tai, tt, utc, tdb };

//! Component order of the displacements of a station loading series.
enum class StationLoadingBasis { unspecified, geocentricRadialEastNorth, cartesian };

//! Sources of loading that a series includes.
enum class StationLoadingComponent { atmosphere, hydrology, ocean };

//! Whether a signal is included in a station loading series.
enum class StationLoadingSignalPresence { unknown, absent, present };

//! Origin of the reference frame of a station loading series.
enum class StationLoadingOrigin { unspecified, wholeEarthCentreOfMass, solidEarthCentreOfMass, centreOfFigure, custom };

//! Epoch of a station loading sample, as given in the file: modified Julian day and seconds of that day.
struct StationLoadingEpoch {
    //! Modified Julian day and seconds since the start of that day.
    int modifiedJulianDay = 0;
    long double secondsOfDay = 0.0L;
};

//! Description of a station loading series. This is supplied by the user; it is not derived from the file name or labels.
struct StationLoadingMetadata {
    //! Name of the file (or other source) and its header lines.
    std::string source;
    std::vector< std::string > sourceHeaders;

    //! Name of the product, its release, and the description of the product on which these settings are based.
    std::string product;
    std::string release;
    std::string sourceProfile;

    //! Reference frame and origin of the displacements, with a text description of the origin.
    std::string referenceFrame;
    StationLoadingOrigin referenceOrigin = StationLoadingOrigin::unspecified;
    std::string originDefinition;

    //! Model of the response of the ocean to atmospheric pressure (for instance an inverted barometer).
    std::string oceanResponse;

    //! Sources of loading included in the series.
    std::vector< StationLoadingComponent > components;

    //! Whether the diurnal, semidiurnal, annual and semiannual signals, the mean and the trend are included in the series.
    StationLoadingSignalPresence diurnal = StationLoadingSignalPresence::unknown;
    StationLoadingSignalPresence semidiurnal = StationLoadingSignalPresence::unknown;
    StationLoadingSignalPresence annual = StationLoadingSignalPresence::unknown;
    StationLoadingSignalPresence semiannual = StationLoadingSignalPresence::unknown;
    StationLoadingSignalPresence mean = StationLoadingSignalPresence::unknown;
    StationLoadingSignalPresence trend = StationLoadingSignalPresence::unknown;

    //! Names of the tidal constituents removed from the series.
    std::vector< std::string > removedHarmonics;

    //! Units of the displacements in the file.
    std::string originalDisplacementUnits = "m";
};

//! Displacement of one site at one epoch of a station loading series.
struct StationLoadingSample {
    //! Displacement (metres), with the components in the order given by the series' basis.
    std::array< double, 3 > displacement{};

    //! Line number and text of the line in the file from which this sample was read.
    std::size_t sourceLine = 0;
    std::string sourceRecord;

    //! Epoch text from columns 10 to 43 of an EPHEDISP line; kept for information only, not used as the sample epoch.
    std::string informationalEpoch;
};

//! Site of a station loading series, with its position and displacement samples.
struct StationLoadingSite {
    //! Name of the site in the file; it is not automatically matched to a ground station name.
    std::string identifier;

    //! Cartesian position of the site (metres), and the distance (metres) from it within which a station may use its data.
    std::array< double, 3 > cartesianPosition{};
    double validityRadiusMetres = 0.0;

    //! Line number and text of the line in the file that defines the site, and the coordinate text of that line.
    std::size_t sourceLine = 0;
    std::string sourceRecord;
    std::string informationalCoordinates;

    //! Samples, with the zero-based index of their epoch in StationLoadingData::epochs as key. Samples may be missing at
    //! the start or end of the series, but not in between.
    std::map< std::size_t, StationLoadingSample > samples;
};

//! Station loading series: displacements of a set of sites at a set of epochs. The epochs may be regular (as in an
//! EPHEDISP file) or irregular. Within each site, the epochs with samples must be consecutive; no gap is filled.
struct StationLoadingData {
    //! Description of the series.
    StationLoadingMetadata metadata;

    //! Time scale of the epochs, and component order of the displacements.
    StationLoadingTimeScale timeScale = StationLoadingTimeScale::tai;
    StationLoadingBasis basis = StationLoadingBasis::unspecified;

    //! Epochs of the samples, in increasing order.
    std::vector< StationLoadingEpoch > epochs;

    //! First and last epoch, and the interval between epochs, as declared in the file header (if given).
    std::optional< StationLoadingEpoch > declaredBeginEpoch;
    std::optional< StationLoadingEpoch > declaredEndEpoch;
    std::optional< long double > nominalSampleIntervalSeconds;

    //! Sites, with their identifier as key.
    std::map< std::string, StationLoadingSite > sites;

    //! Text of the first and last epoch in the file header; kept for information only.
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
