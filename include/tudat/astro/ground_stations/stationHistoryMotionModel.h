/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_STATION_HISTORY_MOTION_MODEL_H
#define TUDAT_STATION_HISTORY_MOTION_MODEL_H
#include "tudat/astro/ground_stations/stationCoordinateHistory.h"
#include "tudat/astro/ground_stations/stationReferencePoint.h"
#include "tudat/astro/ground_stations/groundStationState.h"
#include "tudat/astro/ground_stations/earthDeformationTime.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"

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
     *  \param converter Time scale converter used to convert TDB to TAI (if nullptr, a converter without EOP data is used)
     */
    StationHistoryMotionModel( std::shared_ptr< StationCoordinateHistory > history,
                               std::vector< std::shared_ptr< StationReferencePointOffsetHistory > > referencePointOffsets = {},
                               std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter = nullptr ):
        history_( std::move( history ) ), referencePointOffsets_( std::move( referencePointOffsets ) ), converter_( std::move( converter ) )
    {
        if( !history_ )
        {
            throw std::runtime_error( "Station history motion requires a station coordinate history" );
        }
        validateStationReferencePointOffsets( history_->getMetadata( ), referencePointOffsets_ );
        if( !converter_ )
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
        const Time tai = converter_->getTaiTime( basic_astrodynamics::tdb_scale, time );
        Eigen::Vector6d state = history_->getState( tai );
        for( const auto& offsets : referencePointOffsets_ )
        {
            state += offsets->getStateCorrection( tai );
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
        const Time tai = converter_->getTaiTime( basic_astrodynamics::tdb_scale, time );
        Eigen::Vector3d position = history_->getPosition( tai );
        for( const auto& offsets : referencePointOffsets_ )
        {
            position += offsets->getPositionCorrection( tai );
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
    //! Station coordinate history from which the position is computed.
    std::shared_ptr< StationCoordinateHistory > history_;

    //! Reference point offset histories added to the position of the history, in order.
    std::vector< std::shared_ptr< StationReferencePointOffsetHistory > > referencePointOffsets_;

    //! Time scale converter used to convert TDB to TAI.
    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter_;
};
}  // namespace ground_stations
}  // namespace tudat
#endif
