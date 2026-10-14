/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_TRANSFORMED_STATION_HISTORY_TRAJECTORY_H
#define TUDAT_TRANSFORMED_STATION_HISTORY_TRAJECTORY_H
#include "tudat/astro/ground_stations/stationHistoryTrajectory.h"
#include "tudat/astro/reference_frames/terrestrialFrameTransformation.h"

namespace tudat
{
namespace ground_stations
{
//! Calendar decimal year/rate from the very same source coordinates as the station.
reference_frames::TerrestrialFrameEvaluationEpoch getStationHistoryTerrestrialFrameEpoch( const StationHistoryTime& time );

//! A complete affine state transform; source coefficients and metadata remain inspectable.
/*!
 * The inherited solution/PSD/seasonal records retain source coordinates and provenance.
 * getMetadata describes the evaluated output; getSourceTrajectory retains its input
 * identity. Origin and tide conventions must already match the explicit target.
 * Release alignment is an explicit caller declaration, never a new Helmert datum.
 */
class TransformedStationHistoryTrajectory : public StationHistoryTrajectory
{
public:
    TransformedStationHistoryTrajectory( std::shared_ptr< StationHistoryTrajectory > source,
                                         const input_output::TerrestrialFrameTransformationData& catalogue,
                                         input_output::TerrestrialFrameIdentity target,
                                         std::string targetOrigin,
                                         std::string targetTideSystem,
                                         bool acknowledgeAlignedReleases = false );
    Eigen::Vector3d getPosition( const StationHistoryTime& time ) const override;
    Eigen::Vector6d getState( const StationHistoryTime& time ) const override;
    const StationHistoryMetadata& getMetadata( ) const override
    {
        return outputMetadata_;
    }
    bool hasFittedSeasonalTerms( ) const override
    {
        return source_->hasFittedSeasonalTerms( );
    }
    const std::shared_ptr< StationHistoryTrajectory >& getSourceTrajectory( ) const
    {
        return source_;
    }
    const input_output::TerrestrialFrameTransformationData& getCatalogue( ) const
    {
        return model_.getData( );
    }
    bool getAcknowledgedReleaseAlignment( ) const
    {
        return acknowledgedReleaseAlignment_;
    }

private:
    std::shared_ptr< StationHistoryTrajectory > source_;
    reference_frames::TerrestrialFrameTransformationModel model_;
    StationHistoryMetadata outputMetadata_;
    bool acknowledgedReleaseAlignment_;
};

std::shared_ptr< TransformedStationHistoryTrajectory > createTransformedStationHistoryTrajectory(
        std::shared_ptr< StationHistoryTrajectory > source,
        const input_output::TerrestrialFrameTransformationData& catalogue,
        input_output::TerrestrialFrameIdentity target,
        std::string targetOrigin,
        std::string targetTideSystem,
        bool acknowledgeAlignedReleases = false );
}  // namespace ground_stations
}  // namespace tudat
#endif
