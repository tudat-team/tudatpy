/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_STATION_HISTORY_MOTION_MODEL_H
#define TUDAT_STATION_HISTORY_MOTION_MODEL_H
#include "tudat/astro/ground_stations/stationCoordinateHistory.h"
#include "tudat/astro/ground_stations/stationReferencePoint.h"
#include "tudat/astro/ground_stations/groundStationState.h"
#include "tudat/astro/ground_stations/earthDeformationTime.h"

namespace tudat
{
namespace ground_stations
{

//! Station motion from a station coordinate history (for instance the ITRF positions, velocities, discontinuities,
//! post-seismic and seasonal terms of a station), optionally moved to another point of the station by reference point
//! offsets (for instance from the marker to the antenna reference point). The motion is the resulting position minus the
//! nominal station position.
class StationHistoryMotionModel : public StationMotionModel
{
public:
    //! Constructor
    /*!
     *  Constructor
     *  \param history Station coordinate history
     *  \param referencePointOffsets Offset histories added to the position of the history, in this order (see
     *  validateStationReferencePointOffsets)
     *  \param timeFunction Function converting a TDB epoch to the time used by the history (with rates, if the boolean is
     *  true). If empty, the IGN convention is used, with the time scale converter below.
     *  \param converter Time scale converter used by the IGN convention (if nullptr, a converter without EOP data is used)
     */
    StationHistoryMotionModel( std::shared_ptr< StationCoordinateHistory > history,
                               std::vector< std::shared_ptr< StationReferencePointOffsetHistory > > referencePointOffsets = {},
                               std::function< StationHistoryTime( const Time&, bool ) > timeFunction = nullptr,
                               std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter = nullptr ):
        history_( std::move( history ) ), referencePointOffsets_( std::move( referencePointOffsets ) ),
        timeFunction_( std::move( timeFunction ) ), converter_( std::move( converter ) )
    {
        if( !history_ )
        {
            throw std::runtime_error( "Station history motion requires a station coordinate history" );
        }
        validateStationReferencePointOffsets( history_->getMetadata( ), referencePointOffsets_ );
        if( !converter_ && !timeFunction_ )
        {
            converter_ = basic_astrodynamics::getGeocentricTimeConverter( );
        }
    }

    //! Position and velocity (per TDB second) of the station point, minus the nominal station position.
    Eigen::Vector6d getBodyFixedStationMotion( const Time& time,
                                               const std::shared_ptr< GroundStationState > station,
                                               const std::string& origin ) override
    {
        if( !station )
        {
            throw std::runtime_error( "Station history motion requires a nominal station state" );
        }
        const StationHistoryTime historyTime = sourceTime( time, true );
        Eigen::Vector6d state = history_->getState( historyTime );
        for( const auto& offsets : referencePointOffsets_ )
        {
            state += offsets->getStateCorrection( historyTime );
        }
        state.head< 3 >( ) -= station->getNominalCartesianPosition( );
        return state;
    }

    //! Position of the station point, minus the nominal station position, without requiring rates.
    Eigen::Vector3d getBodyFixedStationPositionCorrection( const Time& time,
                                                           const std::shared_ptr< GroundStationState > station,
                                                           const std::string& origin ) override
    {
        if( !station )
        {
            throw std::runtime_error( "Station history motion requires a nominal station state" );
        }
        const StationHistoryTime historyTime = sourceTime( time, false );
        Eigen::Vector3d position = history_->getPosition( historyTime );
        for( const auto& offsets : referencePointOffsets_ )
        {
            position += offsets->getPositionCorrection( historyTime );
        }
        return position - station->getNominalCartesianPosition( );
    }

    //! Station coordinate history.
    const std::shared_ptr< StationCoordinateHistory >& getHistory( ) const
    {
        return history_;
    }

    //! Reference point offset histories, in the order in which they are added.
    const std::vector< std::shared_ptr< StationReferencePointOffsetHistory > >& getReferencePointOffsets( ) const
    {
        return referencePointOffsets_;
    }

private:
    //! Source coordinates of the TDB epoch, from the custom callback or the IGN nominal UTC profile.
    StationHistoryTime sourceTime( const Time& time, bool rates ) const
    {
        return timeFunction_ ? timeFunction_( time, rates ) : getIgnNominalUtcStationTimeFromTdb( time, converter_, rates );
    }

    //! Station coordinate history from which the position is computed.
    std::shared_ptr< StationCoordinateHistory > history_;

    //! Reference point offset histories added to the position of the history, in order.
    std::vector< std::shared_ptr< StationReferencePointOffsetHistory > > referencePointOffsets_;

    //! Function converting a TDB epoch to the time used by the station history (with rates, if the boolean is true). If
    //! empty, the IGN convention is used, with the time scale converter below.
    std::function< StationHistoryTime( const Time&, bool ) > timeFunction_;
    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter_;
};
}  // namespace ground_stations
}  // namespace tudat
#endif
