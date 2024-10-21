#ifndef TUDAT_POLETIDEDEFORMATION_H
#define TUDAT_POLETIDEDEFORMATION_H


#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/earth_orientation/polarMotionCalculator.h"
#include "tudat/astro/ground_stations/bodyDeformationModel.h"
#include "tudat/interface/sofa/fundamentalArguments.h"


namespace tudat
{

namespace basic_astrodynamics
{

//! This class is used to calculate the displacement of a site on Earth due to pole tide deformation.
/*!
 *  This class is used to calculate the displacement of a site on Earth due to pole tide deformation,
 *  as described in Section 7.1.4 of the IERS 2010 Conventions.
 */
class PoleTideDeformation: public BodyDeformationModel
{
public:

    //! Constructor creating pole tide deformation object from polar motion model.
    /*!
     *  Constructor creating pole tide deformation object from polar motion model.
     *  \param polarMotionCalculator Object capable of calculating pole position in ITRS at a given time instant.
     */
    PoleTideDeformation( const std::shared_ptr< earth_orientation::PolarMotionCalculator > polarMotionCalculator ):
        polarMotionCalculator_( polarMotionCalculator ){ }

    ~PoleTideDeformation( ){ }


    Eigen::Vector3d calculateEnuDisplacement(
        const double ephemerisTime,
        const std::shared_ptr< ground_stations::GroundStationState > nominalSiteState );

    Eigen::Vector3d calculateDisplacement(
        const double ephemerisTime,
        const std::shared_ptr< ground_stations::GroundStationState > stationState );

    Eigen::Vector3d calculateDisplacement(
        const double time,
        const Eigen::Vector3d& bodyFixedPosition )
    {
        return calculateDisplacement( time, std::make_shared< ground_stations::GroundStationState >( bodyFixedPosition, coordinate_conversions::cartesian_position ) );
    }


private:

    //! Pointer to object to calculate the polar motion at a given time.
    /*!
     *  Pointer to object to calculate the polar motion at a given time.
     */
    std::shared_ptr< earth_orientation::PolarMotionCalculator > polarMotionCalculator_;
};


}

}
#endif // TUDAT_POLETIDEDEFORMATION_H
