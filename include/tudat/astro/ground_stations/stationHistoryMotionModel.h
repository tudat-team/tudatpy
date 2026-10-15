/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_STATION_HISTORY_MOTION_MODEL_H
#define TUDAT_STATION_HISTORY_MOTION_MODEL_H
#include "tudat/astro/ground_stations/stationHistoryTrajectory.h"
#include "tudat/astro/ground_stations/groundStationState.h"
#include "tudat/astro/ground_stations/earthDeformationTime.h"

namespace tudat
{
namespace ground_stations
{
class StationHistoryMotionModel : public StationMotionModel
{
public:
    //! The boolean requests rates; position callbacks may leave both rates unavailable.
    StationHistoryMotionModel( std::shared_ptr< StationHistoryTrajectory > trajectory,
                               std::function< StationHistoryTime( const Time&, bool ) > timeFunction = nullptr,
                               std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter = nullptr ):
        trajectory_( std::move( trajectory ) ), timeFunction_( std::move( timeFunction ) ), converter_( std::move( converter ) )
    {
        if( !trajectory_ ) throw std::runtime_error( "Station history motion requires a trajectory" );
        if( !converter_ && !timeFunction_ ) converter_ = basic_astrodynamics::getGeocentricTimeConverter( );
    }

    //! Absolute trajectory state minus the nominal station position; velocity per TDB second.
    Eigen::Vector6d getBodyFixedStationMotion( const Time& time,
                                               const std::shared_ptr< GroundStationState > station,
                                               const std::string& origin ) override
    {
        if( !station )
        {
            throw std::runtime_error( "Station history motion requires a nominal station state" );
        }
        Eigen::Vector6d state = trajectory_->getState( sourceTime( time, true ) );
        state.head< 3 >( ) -= station->getNominalCartesianPosition( );
        return state;
    }

    //! Absolute trajectory position minus the nominal station position, without requiring rates.
    Eigen::Vector3d getBodyFixedStationPositionCorrection( const Time& time,
                                                           const std::shared_ptr< GroundStationState > station,
                                                           const std::string& origin ) override
    {
        if( !station )
        {
            throw std::runtime_error( "Station history motion requires a nominal station state" );
        }
        return trajectory_->getPosition( sourceTime( time, false ) ) - station->getNominalCartesianPosition( );
    }

    const std::shared_ptr< StationHistoryTrajectory >& getTrajectory( ) const
    {
        return trajectory_;
    }

private:
    //! Source coordinates of the TDB epoch, from the custom callback or the IGN nominal UTC profile.
    StationHistoryTime sourceTime( const Time& time, bool rates ) const
    {
        return timeFunction_ ? timeFunction_( time, rates ) : getIgnNominalUtcStationTimeFromTdb( time, converter_, rates );
    }

    std::shared_ptr< StationHistoryTrajectory > trajectory_;
    std::function< StationHistoryTime( const Time&, bool ) > timeFunction_;
    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter_;
};
}  // namespace ground_stations
}  // namespace tudat
#endif
