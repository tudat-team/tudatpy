#include <boost/assign/list_of.hpp>

#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalBodyShapeModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/oblateSpheroidBodyShapeModel.h"
#include "Tudat/Astrodynamics/GroundStations/groundStationState.h"

namespace tudat
{

namespace ground_stations
{

//! Constructor
GroundStationState::GroundStationState(
        const Eigen::Vector3d stationPosition,
        const coordinate_conversions::PositionElementTypes inputElementType,
        const boost::shared_ptr< basic_astrodynamics::BodyShapeModel > bodySurface ):
    bodySurface_( bodySurface )
{
    resetGroundStationPositionAtEpoch( stationPosition, inputElementType );
}

//! Function to obtain the Cartesian state of the ground station in the local frame at a given time.
basic_mathematics::Vector6d GroundStationState::getCartesianStateInTime(
        const double secondsSinceEpoch,
        const double inputReferenceEpoch )
{
    return ( basic_mathematics::Vector6d( ) << cartesianPosition_, Eigen::Vector3d::Zero( ) ).finished( );
}

//! Function to (re)set the nominal state of the station
void GroundStationState::resetGroundStationPositionAtEpoch(
                const Eigen::Vector3d stationPosition,
                const coordinate_conversions::PositionElementTypes inputElementType )
{
    using namespace coordinate_conversions;
    using mathematical_constants::PI;

    // Set Cartesian and spherical position
    cartesianPosition_ = coordinate_conversions::convertPositionElements(
                stationPosition, inputElementType, coordinate_conversions::cartesian_position, bodySurface_ );
    sphericalPosition_ = coordinate_conversions::convertPositionElements(
                stationPosition, inputElementType, coordinate_conversions::spherical_position, bodySurface_ );

    // If possible, set geodetic position, otherwise, set to NaN.
    try
    {
        geodeticPosition = coordinate_conversions::convertPositionElements(
                    stationPosition, inputElementType, coordinate_conversions::geodetic_position, bodySurface_ );
    }
    catch( std::runtime_error )
    {
        geodeticPosition = Eigen::Vector3d::Constant( TUDAT_NAN );
    }

}

}


}
