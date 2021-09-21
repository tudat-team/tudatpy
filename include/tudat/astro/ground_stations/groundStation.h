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
                   const std::string& stationId ):
        nominalStationState_( stationState ),  pointingAnglesCalculator_( pointingAnglesCalculator ), stationId_( stationId ){ }


    //! Function that returns (at reference epoch) the state of the ground station
    /*!
     *  Function that returns (at reference epoch) the state of the ground station. State is computed by nominalStationState_
     *  and cast to the required format by this function
     *  \param time Time (in seconds since J2000) at which state is to be retrieved
     *  \return State at requested time.
     */
    template< typename StateScalarType, typename TimeType >
    Eigen::Matrix< StateScalarType, 6, 1 > getStateInPlanetFixedFrame( const TimeType& time )
    {
        return ( nominalStationState_->getCartesianStateInTime(
                     static_cast< double >( time ), basic_astrodynamics::JULIAN_DAY_ON_J2000  ) ).
                template cast< StateScalarType >( );
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

private:

    //! Object to define and compute the state of the ground station.
    std::shared_ptr< GroundStationState > nominalStationState_;

    //! Object used to computed pointing angles (elevation, azimuth) to a given target from this ground station.
    std::shared_ptr< PointingAnglesCalculator > pointingAnglesCalculator_;

    //! Name of the ground station
    std::string stationId_;
};

} // namespace ground_stations

} // namespace tudat

#endif // TUDAT_GROUNDSTATION_H
