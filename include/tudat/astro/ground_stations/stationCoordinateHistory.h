/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_STATION_COORDINATE_HISTORY_H
#define TUDAT_STATION_COORDINATE_HISTORY_H

#include <functional>
#include <memory>
#include <optional>
#include <Eigen/Geometry>
#include "tudat/basics/timeType.h"
#include "tudat/io/readSinexStationHistory.h"
#include "tudat/io/readTerrestrialFrameTransformationData.h"

namespace tudat
{
namespace earth_orientation
{
class TerrestrialTimeScaleConverter;
}
namespace ground_stations
{
//! Independent source coordinates; physical validity is always selected on TAI.
/*!
 * MJD and calendar year are split to avoid subtracting large floating-point epochs.
 * Rates, when supplied, are derivatives with respect to environment TDB seconds.
 * The coordinate convention must agree with that used to construct the station coordinate history.
 */
struct StationHistoryTime {
    //! The epoch as TAI (seconds since J2000).
    Time tai;

    //! The epoch as modified Julian day and fraction of that day.
    int mjdDay = 0;
    long double mjdFraction = 0.0L;

    //! The epoch as calendar year and fraction of that year.
    int calendarYear = 0;
    long double yearFraction = 0.0L;

    //! Derivatives of the modified Julian date and of the decimal year with respect to TDB, if computed.
    std::optional< double > mjdRate;
    std::optional< double > decimalYearRate;

    //! Name of the convention used to compute the modified Julian date and year from the epoch.
    std::string convention;
};

//! Strict nominal civil UTC profile, post-1972; inserted seconds are not representable.
StationHistoryTime getIgnNominalUtcStationTimeFromTdb(
        const Time& tdb,
        std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter = nullptr,
        bool requireRates = false );

//! Explicit TAI input; a missing clock rate still permits position evaluation.
StationHistoryTime getIgnNominalUtcStationTimeFromTai( const Time& tai, std::optional< double > taiRateWithRespectToTdb = std::nullopt );
StationHistoryTime getIgnNominalUtcStationTime( const input_output::SinexCivilUtc& civil );

//! Physical interval conversion also accepts a real civil inserted-second label.
Time convertSinexCivilUtcToTai( const input_output::SinexCivilUtc& civil );

//! Validity interval of a station history record, in TAI.
struct StationHistoryInterval {
    //! Explicit unbounded ends; finite intervals are [startTai, endTai).
    std::optional< Time > startTai;
    std::optional< Time > endTai;
};

//! Station position and velocity, valid in one time interval (one solution of the station in the source file).
struct StationHistoryLinearSolution {
    //! Solution as read from the source file, and the interval in which its position is valid.
    input_output::SinexCartesianSolution source;
    StationHistoryInterval validity;

    //! Independent V history identifier; SSC velocities remain attached to their P solution.
    std::string velocitySolutionId;
    StationHistoryInterval velocityValidity;

    //! Reference epoch, and the position (metres) and velocity (metres per Julian year) at that epoch.
    StationHistoryTime referenceTime;
    Eigen::Vector3d positionMetres = Eigen::Vector3d::Zero( );
    Eigen::Vector3d velocityMetresPerJulianYear = Eigen::Vector3d::Zero( );
};

//! Function of time of a post-seismic deformation term.
enum class StationPostSeismicFunction { logarithmic, exponential };

//! Post-seismic deformation of a station after an earthquake: a logarithmic or exponential decay with time.
struct StationPostSeismicTerm {
    //! Logarithmic or exponential function, and the epoch of the earthquake.
    StationPostSeismicFunction function = StationPostSeismicFunction::logarithmic;
    StationHistoryTime eventTime;

    //! Amplitude (metres, per Cartesian component) and relaxation time (Julian years) of the term.
    Eigen::Vector3d amplitudeMetres = Eigen::Vector3d::Zero( );
    double relaxationJulianYears = 0.0;

    //! Parameter numbers of the amplitude and relaxation time in the source file, and the line from which the term was read.
    std::size_t amplitudeParameterIndex = 0;
    std::size_t relaxationParameterIndex = 0;
    input_output::SinexSourceRecord provenance;
};

//! Periodic (seasonal) motion of a station, with a given number of cycles per calendar year.
struct StationSeasonalTerm {
    //! Coefficients belong to a position solution, never to a velocity interval ID.
    std::string solutionId;

    //! Number of cycles per calendar year, and the cosine and sine amplitudes (metres, per Cartesian component).
    int cyclesPerCalendarYear = 1;
    Eigen::Vector3d cosineMetres = Eigen::Vector3d::Zero( );
    Eigen::Vector3d sineMetres = Eigen::Vector3d::Zero( );

    //! Parameter numbers of the amplitudes in the source file.
    std::vector< std::size_t > parameterIndices;
};

//! Point of a station whose position is described: the monument marker, the antenna reference point, or a calibrated
//! point of an instrument.
enum class StationReferencePointKind { marker, antennaReferencePoint, calibratedInstrumentPoint };

//! Description of a station history: which station and point it describes, and its reference frame and conventions.
struct StationHistoryMetadata {
    //! Station identifier and DOMES number.
    input_output::SinexSiteKey site;
    std::string domesId;

    //! Reference frame and its release, origin and tide system of the positions.
    std::string frame;
    std::string release;
    std::string origin;
    std::string tideSystem;

    //! Convention used to convert epochs to the time used by the history (see StationHistoryTime).
    std::string timeConvention;

    //! Point of the station whose position is described, and its name.
    StationReferencePointKind referencePointKind = StationReferencePointKind::marker;
    std::string referencePointName = "marker";

    //! Source files of the history, and the kinds of station motion it includes (for instance post-seismic deformation).
    std::vector< std::string > sources;
    std::vector< std::string > signals;
};

//! Catalogued coordinates of a station point over time: linear (position and velocity) solutions with consecutive
//! validity intervals, plus post-seismic and seasonal terms.
/*!
 * All applicable post-seismic terms accumulate. At a genuine event the correction
 * is zero and its rate is the right-hand rate. At a solution boundary the new
 * half-open solution supplies both position and the right-hand velocity. Neither
 * value asserts a two-sided derivative at a discontinuity. Position-only calls
 * never require source-coordinate rates. Gaps, overlaps and extrapolation fail.
 */
class StationCoordinateHistory
{
public:
    //! Own and validate the published solution and correction records once.
    StationCoordinateHistory( std::vector< StationHistoryLinearSolution > solutions,
                              std::vector< StationPostSeismicTerm > postSeismicTerms,
                              std::vector< StationSeasonalTerm > seasonalTerms,
                              StationHistoryMetadata metadata );

    //! Evaluate absolute terrestrial position without requiring source-coordinate rates.
    Eigen::Vector3d getPosition( const StationHistoryTime& time ) const;

    //! Evaluate absolute terrestrial position and its analytic TDB velocity.
    Eigen::Matrix< double, 6, 1 > getState( const StationHistoryTime& time ) const;

    //! Description of the history: station, point, frame and conventions.
    const StationHistoryMetadata& getMetadata( ) const
    {
        return metadata_;
    }

    //! Linear solutions.
    const std::vector< StationHistoryLinearSolution >& getSolutions( ) const
    {
        return solutions_;
    }

    //! Post-seismic terms.
    const std::vector< StationPostSeismicTerm >& getPostSeismicTerms( ) const
    {
        return postSeismicTerms_;
    }

    //! Seasonal terms.
    const std::vector< StationSeasonalTerm >& getSeasonalTerms( ) const
    {
        return seasonalTerms_;
    }

    //! Whether the history contains fitted seasonal terms.
    bool hasFittedSeasonalTerms( ) const
    {
        return !seasonalTerms_.empty( );
    }

private:
    //! Evaluate the linear, cumulative post-seismic and fitted seasonal records.
    Eigen::Matrix< double, 6, 1 > evaluate( const StationHistoryTime& time, bool rates ) const;

    //! Linear solutions, ordered in time.
    std::vector< StationHistoryLinearSolution > solutions_;

    //! Post-seismic terms (all terms whose earthquake precedes an epoch are added) and seasonal terms.
    std::vector< StationPostSeismicTerm > postSeismicTerms_;
    std::vector< StationSeasonalTerm > seasonalTerms_;

    //! Description of the history.
    StationHistoryMetadata metadata_;
};

//! Transform a station coordinate history to another terrestrial reference frame realization.
/*!
 * The linearized Helmert transformation (translation T, rotation and scale matrix M, and their rates) of the catalogue
 * is applied to the coefficients once, and a new history is returned. For each linear solution, with reference epoch t0,
 * position x0 and velocity v, the new position is T(t0) + M(t0) x0 and the new velocity Tdot + Mdot x0 + M(t0) v. The
 * post-seismic and seasonal amplitudes are multiplied by M. Compared with evaluating the transformation at each epoch,
 * this treats the calendar decimal years of the parameter rates as Julian years (they differ by up to about 0.002 year,
 * so the position differs by up to a few micrometres for the oldest ITRF realizations), and omits the term
 * Mdot v (t - t0)^2 (tens of nanometres).
 * The reference origin and tide system are not converted, and must be equal to those of the input history.
 * \param history Station coordinate history that is to be transformed
 * \param catalogue Table of transformation parameters
 * \param target Realization and release of the target frame
 * \param targetOrigin Reference origin of the output (must be equal to that of the input)
 * \param targetTideSystem Tide system of the output (must be equal to that of the input)
 * \param acknowledgeAlignedReleases Whether a release that differs from that of the catalogue (for instance an update of
 * a realization) is accepted as aligned with it. This does not add a transformation.
 * \return Station coordinate history in the target frame
 */
std::shared_ptr< StationCoordinateHistory > transformStationCoordinateHistory(
        const std::shared_ptr< StationCoordinateHistory >& history,
        const input_output::TerrestrialFrameTransformationData& catalogue,
        input_output::TerrestrialFrameIdentity target,
        const std::string& targetOrigin,
        const std::string& targetTideSystem,
        bool acknowledgeAlignedReleases = false );

//! Explicit profile for IGN's source-index pairing and relaxation-unit convention.
/*!
 * Every product must carry consistent supplied frame/release/origin/tide metadata.
 * An ENU product requires an explicit, constant proper ENU-to-terrestrial rotation.
 * The native profile interprets TLOG/TEXP as Julian years even where IGN labels m;
 * the neutral input records remain unchanged. Amplitude index i is paired only
 * with matching relaxation index i+1; input line order has no meaning.
 * SSC velocity coefficients stay attached to their P solution. Independent V
 * intervals provide their own identity and must contain each complete P interval.
 * Shared V intervals require componentwise agreement within the caller's declared
 * tolerance (default 10 nm per Julian year); supplied coefficients are unchanged.
 * A DOMES selector resolves duplicated SITE/ID metadata; the caller warrants the
 * association of parameter records to that monument using companion products.
 */
std::shared_ptr< StationCoordinateHistory > createIgnStationCoordinateHistory(
        const input_output::SinexStationHistory& solutions,
        const input_output::SinexStationHistory& discontinuities,
        const input_output::SinexSiteKey& site,
        const input_output::SinexStationHistory* postSeismic = nullptr,
        const input_output::SinexStationHistory* seasonal = nullptr,
        const std::optional< Eigen::Matrix3d >& enuToTerrestrial = std::nullopt,
        const std::optional< std::string >& domesSelector = std::nullopt,
        double velocityConsistencyToleranceMetresPerJulianYear = 1.0e-8 );

}  // namespace ground_stations
}  // namespace tudat
#endif
