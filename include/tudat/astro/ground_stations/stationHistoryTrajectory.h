/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_STATION_HISTORY_TRAJECTORY_H
#define TUDAT_STATION_HISTORY_TRAJECTORY_H

#include <functional>
#include <memory>
#include <optional>
#include <Eigen/Geometry>
#include "tudat/basics/timeType.h"
#include "tudat/io/readSinexStationHistory.h"

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
 * The coordinate convention must agree with that used to construct the trajectory.
 */
struct StationHistoryTime {
    Time tai;
    int mjdDay = 0;
    long double mjdFraction = 0.0L;
    int calendarYear = 0;
    long double yearFraction = 0.0L;
    std::optional< double > mjdRate;
    std::optional< double > decimalYearRate;
    std::string convention;
};

//! Strict nominal civil UTC profile, post-1972; inserted seconds are not representable.
StationHistoryTime getIgnNominalUtcStationTimeFromTdb(
        const Time& tdb,
        std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter = nullptr,
        bool requireRates = false );
StationHistoryTime getIgnNominalUtcStationTimeFromTdb(
        double tdb,
        std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter = nullptr,
        bool requireRates = false );
//! Explicit TAI input; a missing clock rate still permits position evaluation.
StationHistoryTime getIgnNominalUtcStationTimeFromTai( const Time& tai, std::optional< double > taiRateWithRespectToTdb = std::nullopt );
StationHistoryTime getIgnNominalUtcStationTime( const input_output::SinexCivilUtc& civil );
//! Physical interval conversion also accepts a real civil inserted-second label.
Time convertSinexCivilUtcToTai( const input_output::SinexCivilUtc& civil );

struct StationHistoryInterval {
    //! Explicit unbounded ends; finite intervals are [startTai, endTai).
    std::optional< Time > startTai;
    std::optional< Time > endTai;
};
struct StationHistoryLinearSolution {
    input_output::SinexCartesianSolution source;
    StationHistoryInterval validity;
    //! Independent V history identifier; SSC velocities remain attached to their P solution.
    std::string velocitySolutionId;
    StationHistoryInterval velocityValidity;
    StationHistoryTime referenceTime;
    Eigen::Vector3d positionMetres = Eigen::Vector3d::Zero( );
    Eigen::Vector3d velocityMetresPerJulianYear = Eigen::Vector3d::Zero( );
};
enum class StationPostSeismicFunction { logarithmic, exponential };
struct StationPostSeismicTerm {
    StationPostSeismicFunction function = StationPostSeismicFunction::logarithmic;
    StationHistoryTime eventTime;
    Eigen::Vector3d amplitudeMetres = Eigen::Vector3d::Zero( );
    double relaxationJulianYears = 0.0;
    std::size_t amplitudeParameterIndex = 0;
    std::size_t relaxationParameterIndex = 0;
    input_output::SinexSourceRecord provenance;
};
struct StationSeasonalTerm {
    //! Coefficients belong to a position solution, never to a velocity interval ID.
    std::string solutionId;
    int cyclesPerCalendarYear = 1;
    Eigen::Vector3d cosineMetres = Eigen::Vector3d::Zero( );
    Eigen::Vector3d sineMetres = Eigen::Vector3d::Zero( );
    std::vector< std::size_t > parameterIndices;
};
enum class StationReferencePointKind { marker, antennaReferencePoint, calibratedInstrumentPoint };
struct StationHistoryMetadata {
    input_output::SinexSiteKey site;
    std::string domesId;
    std::string frame;
    std::string release;
    std::string origin;
    std::string tideSystem;
    std::string timeConvention;
    StationReferencePointKind referencePointKind = StationReferencePointKind::marker;
    std::string referencePointName = "marker";
    std::vector< std::string > sources;
    std::vector< std::string > signals;
};

//! Absolute terrestrial trajectory with analytic TDB velocity and explicit source time.
/*!
 * All applicable post-seismic terms accumulate. At a genuine event the correction
 * is zero and its rate is the right-hand rate. At a solution boundary the new
 * half-open solution supplies both position and the right-hand velocity. Neither
 * value asserts a two-sided derivative at a discontinuity. Position-only calls
 * never require source-coordinate rates. Gaps, overlaps and extrapolation fail.
 */
class StationHistoryTrajectory
{
public:
    //! Own and validate the published solution and correction records once.
    StationHistoryTrajectory( std::vector< StationHistoryLinearSolution > solutions,
                              std::vector< StationPostSeismicTerm > postSeismicTerms,
                              std::vector< StationSeasonalTerm > seasonalTerms,
                              StationHistoryMetadata metadata );

    //! Permit ordered trajectory wrappers to be destroyed through their base type.
    virtual ~StationHistoryTrajectory( ) = default;

    //! Evaluate absolute terrestrial position without requiring source-coordinate rates.
    virtual Eigen::Vector3d getPosition( const StationHistoryTime& time ) const;

    //! Evaluate absolute terrestrial position and its analytic TDB velocity.
    virtual Eigen::Matrix< double, 6, 1 > getState( const StationHistoryTime& time ) const;

    //! Return the evaluated trajectory's frame, monument and reference-point identity.
    virtual const StationHistoryMetadata& getMetadata( ) const
    {
        return source_ ? source_->getMetadata( ) : metadata_;
    }

    //! Inspect the original published solutions through any ordered wrappers.
    const std::vector< StationHistoryLinearSolution >& getSolutions( ) const
    {
        return source_ ? source_->getSolutions( ) : solutions_;
    }

    //! Inspect the original cumulative post-seismic coefficients without copying them.
    const std::vector< StationPostSeismicTerm >& getPostSeismicTerms( ) const
    {
        return source_ ? source_->getPostSeismicTerms( ) : postSeismicTerms_;
    }

    //! Inspect the original fitted seasonal coefficients without copying them.
    const std::vector< StationSeasonalTerm >& getSeasonalTerms( ) const
    {
        return source_ ? source_->getSeasonalTerms( ) : seasonalTerms_;
    }

    //! Identify fitted seasonal motion for loading-overlap decisions.
    virtual bool hasFittedSeasonalTerms( ) const
    {
        return source_ ? source_->hasFittedSeasonalTerms( ) : !seasonalTerms_.empty( );
    }

protected:
    //! Retain an already validated source for an ordered frame or reference-point wrapper.
    explicit StationHistoryTrajectory( std::shared_ptr< StationHistoryTrajectory > source );

    //! Shared source lifetime; null only for trajectories owning their published records.
    const std::shared_ptr< StationHistoryTrajectory > source_;

private:
    //! Evaluate the owned linear, cumulative post-seismic and fitted seasonal records.
    Eigen::Matrix< double, 6, 1 > evaluate( const StationHistoryTime& time, bool rates ) const;

    //! Published coefficients retained only by the original trajectory.
    std::vector< StationHistoryLinearSolution > solutions_;

    std::vector< StationPostSeismicTerm > postSeismicTerms_;

    std::vector< StationSeasonalTerm > seasonalTerms_;

    //! Source coordinate conventions and provenance for the owned records.
    StationHistoryMetadata metadata_;

};

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
std::shared_ptr< StationHistoryTrajectory > createIgnStationHistoryTrajectory(
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
