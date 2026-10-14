/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_TRANSFORMED_STATION_HISTORY_TRAJECTORY_H
#define TUDAT_TRANSFORMED_STATION_HISTORY_TRAJECTORY_H
#include "tudat/astro/ground_stations/stationHistoryTrajectory.h"
#include "tudat/astro/reference_frames/terrestrialFrameTransformation.h"

namespace tudat
{
namespace ground_stations
{
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
    //! Retain the validated source and bind a catalogue to explicit target coordinate conventions.
    TransformedStationHistoryTrajectory( std::shared_ptr< StationHistoryTrajectory > source,
                                         const input_output::TerrestrialFrameTransformationData& catalogue,
                                         input_output::TerrestrialFrameIdentity target,
                                         std::string targetOrigin,
                                         std::string targetTideSystem,
                                         bool acknowledgeAlignedReleases = false );

    //! Apply the affine position map at the source calendar decimal year.
    Eigen::Vector3d getPosition( const StationHistoryTime& time ) const override;

    //! Include parameter derivatives using the source calendar-year rate with respect to TDB.
    Eigen::Vector6d getState( const StationHistoryTime& time ) const override;

    //! Describe output coordinates while retaining the original monument and reference-point identity.
    const StationHistoryMetadata& getMetadata( ) const override
    {
        return outputMetadata_;
    }

    //! Inspect the retained source and its input-coordinate identity.
    const std::shared_ptr< StationHistoryTrajectory >& getSourceTrajectory( ) const
    {
        return source_;
    }

    //! Inspect the validated transformation parameters and provenance.
    const input_output::TerrestrialFrameTransformationData& getCatalogue( ) const
    {
        return model_.getData( );
    }

    //! Report the caller's acknowledgement of aligned catalogue releases.
    bool getAcknowledgedReleaseAlignment( ) const
    {
        return acknowledgedReleaseAlignment_;
    }

private:
    //! Validated affine catalogue, evaluated-output identity and release-alignment declaration.
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
