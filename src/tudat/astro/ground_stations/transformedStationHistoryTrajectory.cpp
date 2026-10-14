/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/astro/ground_stations/transformedStationHistoryTrajectory.h"
#include <cmath>
#include <stdexcept>

namespace tudat
{
namespace ground_stations
{
namespace
{
//! Recover the published calendar year from the trajectory's split coordinate.
double decimalYear( const StationHistoryTime& time )
{
    if( !std::isfinite( time.yearFraction ) || time.yearFraction < 0 || time.yearFraction >= 1 )
    {
        throw std::runtime_error( "Invalid station-history calendar coordinates for frame transformation" );
    }
    return static_cast< double >( time.calendarYear + time.yearFraction );
}

//! Resolve the release identity explicitly supplied by the selected catalogue row.
std::string catalogueRelease( const input_output::TerrestrialFrameTransformationData& catalogue, const std::string& realization )
{
    if( realization == catalogue.referenceFrame.realization )
    {
        return catalogue.referenceFrame.release;
    }

    // The row can declare a realization without a separate release label.
    const auto& identity = catalogue.getParameters( realization ).targetFrame;

    // The published target rows often declare only the realization. An aligned
    // update must be acknowledged explicitly rather than inferred from that omission.
    return identity.release.empty( ) ? identity.realization : identity.release;
}
}  // namespace

//! Validate target coordinate semantics while retaining the source monument and point identity.
TransformedStationHistoryTrajectory::TransformedStationHistoryTrajectory( std::shared_ptr< StationHistoryTrajectory > source,
                                                                          const input_output::TerrestrialFrameTransformationData& catalogue,
                                                                          input_output::TerrestrialFrameIdentity target,
                                                                          std::string targetOrigin,
                                                                          std::string targetTideSystem,
                                                                          bool acknowledgeAlignedReleases ):
    StationHistoryTrajectory( std::move( source ) ),
    model_( catalogue ), outputMetadata_( source_->getMetadata( ) ),
    acknowledgedReleaseAlignment_( acknowledgeAlignedReleases )
{
    // Origin and tide-system declarations must agree before any frame map is applied.
    const auto& input = source_->getMetadata( );

    if( target.realization.empty( ) || target.release.empty( ) || input.origin.empty( ) || input.tideSystem.empty( ) ||
        targetOrigin != input.origin || targetTideSystem != input.tideSystem )
    {
        throw std::runtime_error( "Station frame transformation requires explicit matching origin/tide conventions and target identity" );
    }
    if( input.frame == target.realization )
    {
        if( input.release != target.release && !acknowledgeAlignedReleases )
        {
            throw std::runtime_error( "Same-frame station release change requires an explicit alignment declaration" );
        }
    }
    else
    {
        // Distinct realizations must use the published release identities or an explicit alignment declaration.
        const auto inputRelease = catalogueRelease( catalogue, input.frame );

        const auto targetRelease = catalogueRelease( catalogue, target.realization );

        if( ( input.release != inputRelease || target.release != targetRelease ) && !acknowledgeAlignedReleases )
        {
            throw std::runtime_error( "Station releases differ from frame-catalogue identities; explicit alignment declaration required" );
        }
    }
    outputMetadata_.frame = std::move( target.realization );
    outputMetadata_.release = std::move( target.release );
    outputMetadata_.sources.push_back( catalogue.provenance.source );
    outputMetadata_.signals.push_back( input.frame == outputMetadata_.frame ? "declared same-frame identity"
                                                                            : "terrestrial affine frame transformation" );
}

//! Transform available position using the trajectory's calendar year, without requesting a rate.
Eigen::Vector3d TransformedStationHistoryTrajectory::getPosition( const StationHistoryTime& time ) const
{
    // Absolute Cartesian metres and the declared source frame determine the affine position map.
    const Eigen::Vector3d position = source_->getPosition( time );

    const auto& input = source_->getMetadata( );

    if( input.frame == outputMetadata_.frame )
    {
        return position;
    }
    return model_.transformPosition( position, input.frame, outputMetadata_.frame, decimalYear( time ) );
}

//! Include the affine parameter derivatives in metres per TDB second.
Eigen::Vector6d TransformedStationHistoryTrajectory::getState( const StationHistoryTime& time ) const
{
    // Use the same source state and calendar-coordinate mapping as position evaluation.
    const Eigen::Vector6d state = source_->getState( time );

    const auto& input = source_->getMetadata( );

    if( input.frame == outputMetadata_.frame )
    {
        return state;
    }
    if( !time.decimalYearRate )
    {
        throw std::runtime_error( "Station-history frame state requires calendar-year derivative with respect to TDB" );
    }
    return model_.transformState( state, input.frame, outputMetadata_.frame, decimalYear( time ), *time.decimalYearRate );
}

//! Construct the affine wrapper from the caller's explicit target frame and provenance declarations.
std::shared_ptr< TransformedStationHistoryTrajectory > createTransformedStationHistoryTrajectory(
        std::shared_ptr< StationHistoryTrajectory > source,
        const input_output::TerrestrialFrameTransformationData& catalogue,
        input_output::TerrestrialFrameIdentity target,
        std::string targetOrigin,
        std::string targetTideSystem,
        bool acknowledgeAlignedReleases )
{
    return std::make_shared< TransformedStationHistoryTrajectory >( std::move( source ),
                                                                    catalogue,
                                                                    std::move( target ),
                                                                    std::move( targetOrigin ),
                                                                    std::move( targetTideSystem ),
                                                                    acknowledgeAlignedReleases );
}
}  // namespace ground_stations
}  // namespace tudat
