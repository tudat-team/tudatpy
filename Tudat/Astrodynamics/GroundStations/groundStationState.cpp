#include <boost/assign/list_of.hpp>

#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/geodeticCoordinateConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalBodyShapeModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/oblateSpheroidBodyShapeModel.h"
#include "Tudat/Astrodynamics/GroundStations/groundStationState.h"

namespace tudat
{

namespace ground_stations
{

GroundStationState::GroundStationState(
        const Eigen::Vector3d stationPosition,
        const coordinate_conversions::PositionElementTypes inputElementType,
        const boost::shared_ptr< basic_astrodynamics::BodyShapeModel > bodySurface ):
    bodySurface_( bodySurface )
{

    resetGroundStationPositionAtEpoch( stationPosition, inputElementType );
}


Eigen::Vector3d GroundStationState::getCartesianPositionInTime(
        const double secondsSinceEpoch,
        const double inputReferenceEpoch )
{
    return cartesianPosition_;
}

void GroundStationState::resetGroundStationPositionAtEpoch(
                const Eigen::Vector3d stationPosition,
                const coordinate_conversions::PositionElementTypes inputElementType )
{
    using namespace coordinate_conversions;
    using mathematical_constants::PI;

    cartesianPosition_ = coordinate_conversions::convertPositionElements(
                stationPosition, inputElementType, coordinate_conversions::cartesian_position, bodySurface_ );
    sphericalPosition_ = coordinate_conversions::convertPositionElements(
                stationPosition, inputElementType, coordinate_conversions::spherical_position, bodySurface_ );
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
