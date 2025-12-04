/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_GROUNDSTATION_H
#define TUDAT_GROUNDSTATION_H

#include <memory>

#include <Eigen/Core>

#include "tudat/astro/ground_stations/groundStationState.h"
#include "tudat/astro/ground_stations/pointingAnglesCalculator.h"
#include "tudat/astro/ground_stations/meteorologicalConditions.h"
#include "tudat/astro/system_models/timingSystem.h"
#include "tudat/astro/ground_stations/transmittingFrequencies.h"
#include "tudat/astro/system_models/vehicleSystems.h"

namespace tudat
{

namespace ground_stations
{

//! Class to store properties of a ground station (i.e. reference point with associated systems on a celestial body)
class GroundStation
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param stationState Object to define and compute the state of the ground station.
     * \param pointingAnglesCalculator Object used to computed pointing angles (elevation, azimuth) to a given target from this
     * ground station.
     * \param stationId Name of the ground station
     */
    GroundStation( const std::shared_ptr< GroundStationState > stationState,
                   const std::shared_ptr< PointingAnglesCalculator > pointingAnglesCalculator,
                   const std::string& stationId,
                   const std::shared_ptr< StationFrequencyInterpolator > transmittingFrequencyCalculator = nullptr ):
        nominalStationState_( stationState ), pointingAnglesCalculator_( pointingAnglesCalculator ), stationId_( stationId ),
        transmittingFrequencyCalculator_( transmittingFrequencyCalculator )
    {
        stationState->setSiteId( stationId );
    }

    //! Function that returns (at reference epoch) the state of the ground station
    /*!
     *  Function that returns (at reference epoch) the state of the ground station. State is computed by nominalStationState_
     *  and cast to the required format by this function
     *  \param time Time (in seconds since J2000) at which state is to be retrieved
     *  \return State at requested time.
     */
    template< typename StateScalarType, typename TimeType >
    Eigen::Matrix< StateScalarType, 6, 1 > getStateInPlanetFixedFrame( const TimeType& time, const std::string& targetFrameOrigin )
    {
        return ( nominalStationState_->getCartesianStateInTime( static_cast< double >( time ), targetFrameOrigin ) )
                .template cast< StateScalarType >( );
    }

    //! Function to return object to define and compute the state of the ground station.
    /*!
     * Function to return object to define and compute the state of the ground station.
     * \return Object to define and compute the state of the ground station.
     */
    std::shared_ptr< GroundStationState > getNominalStationState( )
    {
        return nominalStationState_;
    }

    //! Function to return name of the ground station
    /*!
     * Function to return name of the ground station
     * \return Name of the ground station
     */
    std::string getStationId( )
    {
        return stationId_;
    }

    //! Function to return object used to computed pointing angles (elevation, azimuth) to a given target from this ground station.
    /*!
     * Function to object used to computed pointing angles (elevation, azimuth) to a given target from this ground station.
     * \return Object used to computed pointing angles (elevation, azimuth) to a given target from this ground station.
     */
    std::shared_ptr< PointingAnglesCalculator > getPointingAnglesCalculator( )
    {
        return pointingAnglesCalculator_;
    }

    std::shared_ptr< system_models::TimingSystem > getTimingSystem( )
    {
        return timingSystem_;
    }

    void setTimingSystem( const std::shared_ptr< system_models::TimingSystem > timingSystem )
    {
        timingSystem_ = timingSystem;
    }

    //! Function to return the object used to compute the ground station's transmitting frequency at a given time
    std::shared_ptr< StationFrequencyInterpolator > getTransmittingFrequencyCalculator( )
    {
        if( transmittingFrequencyCalculator_ == nullptr )
        {
            throw std::runtime_error( "Error when retrieving the frequency calculator for ground station " + stationId_ +
                                      ": no frequency calculator has been defined" );
        }

        return transmittingFrequencyCalculator_;
    }

    bool hasFrequencyCalculator( )
    {
        return ( !( transmittingFrequencyCalculator_ == nullptr  ) ) ;
    }

    //! Function to set the object used to compute the ground station's transmitting frequency at a given time
    void setTransmittingFrequencyCalculator( std::shared_ptr< StationFrequencyInterpolator > transmittingFrequencyCalculator )
    {
        transmittingFrequencyCalculator_ = transmittingFrequencyCalculator;
    }

    std::function< double( const double time ) > getTemperatureFunction( )
    {
        if( meteoData_ == nullptr )
        {
            throw std::runtime_error( "Error when getting temperature function from ground station " + stationId_ +
                                      ": function is not defined." );
        }
        return std::bind( &StationMeteoData::getTemperature, meteoData_, std::placeholders::_1 );
    }

    std::function< double( const double time ) > getPressureFunction( )
    {
        if( meteoData_ == nullptr )
        {
            throw std::runtime_error( "Error when getting pressure function from ground station " + stationId_ +
                                      ": function is not defined." );
        }
        return std::bind( &StationMeteoData::getPressure, meteoData_, std::placeholders::_1 );
    }

    std::function< double( const double time ) > getWaterVaporPartialPressureFunction( )
    {
        if( meteoData_ == nullptr )
        {
            throw std::runtime_error( "Error when getting water vapor partial pressure function from ground station " + stationId_ +
                                      ": function is not defined." );
        }
        return std::bind( &StationMeteoData::getWaterVaporPartialPressure, meteoData_, std::placeholders::_1 );
    }

    std::function< double( const double time ) > getRelativeHumidityFunction( )
    {
        if( meteoData_ == nullptr )
        {
            throw std::runtime_error( "Error when getting relative humidity function from ground station " + stationId_ +
                                      ": function is not defined." );
        }
        return std::bind( &StationMeteoData::getRelativeHumidity, meteoData_, std::placeholders::_1 );
    }

    std::function< double( const double time ) > getDewPointFunction( )
    {
        if( meteoData_ == nullptr )
        {
            throw std::runtime_error( "Error when getting dew point function from ground station " + stationId_ +
                                      ": function is not defined." );
        }
        return std::bind( &StationMeteoData::getDewPointTemperature, meteoData_, std::placeholders::_1 );
    }

    //! Function to retrieve container object with hardware systems present on/in body
    /*!
     * Function to retrieve container object with hardware systems present on/in body.
     * \return Container object with hardware systems present on/in body.
     */
    std::shared_ptr< system_models::VehicleSystems > getVehicleSystems( )
    {
        if( vehicleSystems_ == nullptr )
        {
            vehicleSystems_ = std::make_shared< system_models::VehicleSystems >( );
        }
        return vehicleSystems_;
    }

    //! Function to set container object with hardware systems present on/in body
    /*!
     * Function to set container object with hardware systems present on/in body (typically only non-nullptr for a vehicle).
     * \param vehicleSystems Container object with hardware systems present on/in body.
     */
    void setVehicleSystems( const std::shared_ptr< system_models::VehicleSystems > vehicleSystems )
    {
        vehicleSystems_ = vehicleSystems;
    }

    std::shared_ptr< StationMeteoData > getMeteoData( )
    {
        return meteoData_;
    }

    void setMeteoData( const std::shared_ptr< StationMeteoData > meteoData )
    {
        meteoData_ = meteoData;
    }

    std::shared_ptr< StationTroposphereData > getTroposphereData( )
    {
        return troposphereData_;
    }

    void setTroposphereData( const std::shared_ptr< StationTroposphereData > troposphereData )
    {
        troposphereData_ = troposphereData;
    }

private:
    //! Object to define and compute the state of the ground station.
    std::shared_ptr< GroundStationState > nominalStationState_;

    //! Object used to computed pointing angles (elevation, azimuth) to a given target from this ground station.
    std::shared_ptr< PointingAnglesCalculator > pointingAnglesCalculator_;

    std::shared_ptr< system_models::TimingSystem > timingSystem_;

    //! Name of the ground station
    std::string stationId_;

    //! Object used to defined and compute the ground station's transmitting frequency.
    std::shared_ptr< StationFrequencyInterpolator > transmittingFrequencyCalculator_;

    std::shared_ptr< StationMeteoData > meteoData_;

    std::shared_ptr< StationTroposphereData > troposphereData_;

    //! Container object with hardware systems present on/in body (typically only non-nullptr for a vehicle).
    std::shared_ptr< system_models::VehicleSystems > vehicleSystems_;
};

}  // namespace ground_stations

}  // namespace tudat

#endif  // TUDAT_GROUNDSTATION_H
