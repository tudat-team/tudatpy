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
//! The boolean requests rates; position callbacks may leave both rates unavailable.
using StationHistoryTimeFunction = std::function< StationHistoryTime( const Time&, bool ) >;

class StationHistoryMotionModel : public StationMotionModel
{
public:
    StationHistoryMotionModel( std::shared_ptr< StationHistoryTrajectory > trajectory,
                               StationHistoryTimeFunction timeFunction = nullptr,
                               std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter = nullptr ):
        trajectory_( std::move( trajectory ) ), timeFunction_( std::move( timeFunction ) ), converter_( std::move( converter ) )
    {
        if( !trajectory_ ) throw std::runtime_error( "Station history motion requires a trajectory" );
        if( !converter_ && !timeFunction_ ) converter_ = basic_astrodynamics::getGeocentricTimeConverter( );
    }
    Eigen::Vector6d getBodyFixedStationMotion( const Time& time,
                                               const std::shared_ptr< GroundStationState > station,
                                               const std::string& origin ) override
    {
        return stateCorrection( time, station );
    }
    Eigen::Vector6d getBodyFixedStationMotion( double time,
                                               const std::shared_ptr< GroundStationState > station,
                                               const std::string& origin ) override
    {
        return stateCorrection( time, station );
    }
    Eigen::Vector3d getBodyFixedStationPositionCorrection( const Time& time,
                                                           const std::shared_ptr< GroundStationState > station,
                                                           const std::string& origin ) override
    {
        return positionCorrection( time, station );
    }
    Eigen::Vector3d getBodyFixedStationPositionCorrection( double time,
                                                           const std::shared_ptr< GroundStationState > station,
                                                           const std::string& origin ) override
    {
        return positionCorrection( time, station );
    }
    const std::shared_ptr< StationHistoryTrajectory >& getTrajectory( ) const
    {
        return trajectory_;
    }

    bool hasFittedSeasonalTerms( ) const override
    {
        return trajectory_->hasFittedSeasonalTerms( );
    }

private:
    template< typename Epoch >
    StationHistoryTime sourceTime( const Epoch& time, bool rates ) const
    {
        return timeFunction_ ? timeFunction_( Time( time ), rates ) : getIgnNominalUtcStationTimeFromTdb( time, converter_, rates );
    }
    template< typename Epoch >
    Eigen::Vector6d stateCorrection( const Epoch& time, const std::shared_ptr< GroundStationState >& station ) const
    {
        if( !station ) throw std::runtime_error( "Station history motion requires a nominal station state" );
        Eigen::Vector6d state = trajectory_->getState( sourceTime( time, true ) );
        state.head< 3 >( ) -= station->getNominalCartesianPosition( );
        return state;
    }
    template< typename Epoch >
    Eigen::Vector3d positionCorrection( const Epoch& time, const std::shared_ptr< GroundStationState >& station ) const
    {
        if( !station ) throw std::runtime_error( "Station history motion requires a nominal station state" );
        return trajectory_->getPosition( sourceTime( time, false ) ) - station->getNominalCartesianPosition( );
    }

    std::shared_ptr< StationHistoryTrajectory > trajectory_;
    StationHistoryTimeFunction timeFunction_;
    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter_;
};
}  // namespace ground_stations
}  // namespace tudat
#endif
