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
using StationHistoryTimeFunction = std::function< StationHistoryTime( const basic_astrodynamics::EarthDeformationEpoch&, bool ) >;

class StationHistoryMotionModel : public StationMotionModel
{
public:
    StationHistoryMotionModel( std::shared_ptr< StationHistoryTrajectory > trajectory,
                               StationHistoryTimeFunction timeFunction = nullptr,
                               std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter = nullptr ):
        trajectory_( std::move( trajectory ) ), timeFunction_( std::move( timeFunction ) ), converter_( std::move( converter ) )
    {
        if( !trajectory_ ) throw std::runtime_error( "Station history motion requires a trajectory" );
        if( !timeFunction_ ) timeFunction_ = []( const auto& epoch, bool rates ) { return getIgnNominalUtcStationTime( epoch, rates ); };
    }
    Eigen::Vector6d getBodyFixedStationMotion( double time,
                                               std::shared_ptr< GroundStationState > station,
                                               const std::string& targetFrameOrigin ) override
    {
        if( !station ) throw std::runtime_error( "Station history motion requires a nominal station state" );
        Eigen::Vector6d state =
                trajectory_->getState( timeFunction_( basic_astrodynamics::EarthDeformationEpoch( time, converter_ ), true ) );
        state.head< 3 >( ) -= station->getNominalCartesianPosition( );
        return state;
    }
    Eigen::Vector3d getBodyFixedStationPositionCorrection( double time,
                                                           std::shared_ptr< GroundStationState > station,
                                                           const std::string& targetFrameOrigin ) override
    {
        if( !station ) throw std::runtime_error( "Station history motion requires a nominal station state" );
        return trajectory_->getPosition( timeFunction_( basic_astrodynamics::EarthDeformationEpoch( time, converter_ ), false ) ) -
                station->getNominalCartesianPosition( );
    }
    const std::shared_ptr< StationHistoryTrajectory >& getTrajectory( ) const
    {
        return trajectory_;
    }

private:
    std::shared_ptr< StationHistoryTrajectory > trajectory_;
    StationHistoryTimeFunction timeFunction_;
    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter_;
};
}  // namespace ground_stations
}  // namespace tudat
#endif
