#include <boost/assign/list_of.hpp>

#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/geodeticCoordinateConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalBodyShapeModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/oblateSpheroidBodyShapeModel.h"
#include "Tudat/Astrodynamics/GroundStations/nominalGroundStationState.h"

namespace tudat
{

namespace ground_stations
{

NominalGroundStationState::NominalGroundStationState(
        const Eigen::Vector3d stationCartesianPosition,
        const boost::shared_ptr< basic_astrodynamics::BodyShapeModel > bodySurface,
        const double referenceJulianYear, const bool setTransformation ):
    bodySurface_( bodySurface ),
    positionVariationsFunction_( boost::lambda::constant( Eigen::Vector3d::Zero( ) ) ),
    referenceJulianYear_( referenceJulianYear )
{
    using namespace coordinate_conversions;
    using mathematical_constants::PI;

    cartesianPosition_ = stationCartesianPosition;
    sphericalPosition_ = convertCartesianToSpherical( cartesianPosition_ );

}


Eigen::Vector3d NominalGroundStationState::getCartesianPositionInTime(
        const double secondsSinceEpoch,
        const double inputReferenceEpoch )
{
    return cartesianPosition_;
}

void NominalGroundStationState::resetGroundStationPositionAtEpoch( const Eigen::Vector3d cartesianPosition )
{
    using namespace coordinate_conversions;
    using mathematical_constants::PI;

    cartesianPosition_ = cartesianPosition;
    sphericalPosition_ = convertCartesianToSpherical( cartesianPosition_ );

}

}


}
