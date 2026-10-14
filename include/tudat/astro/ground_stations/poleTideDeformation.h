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
 *  as described in Section 7.1.4 of the IERS Conventions, including the secular pole
 *  update of 1 February 2018. Pole tide amplitudes use geocentric station latitude.
 */
class PoleTideDeformation : public BodyDeformationModel
{
public:
    using BodyDeformationModel::calculateDisplacementState;
    using BodyDeformationModel::calculateDisplacementVelocity;
    Eigen::Vector3d calculateDisplacement( const EarthDeformationEpoch& epoch, const Eigen::Vector3d& position ) override;
    Eigen::Vector3d calculateDisplacement( const EarthDeformationEpoch& epoch,
                                           const std::shared_ptr< ground_stations::GroundStationState > station ) override;
    Eigen::Vector3d calculateDisplacementVelocity( const EarthDeformationEpoch& epoch, const Eigen::Vector3d& position ) override;
    Eigen::Vector3d calculateDisplacementVelocity( const EarthDeformationEpoch& epoch,
                                                   const std::shared_ptr< ground_stations::GroundStationState > station ) override;
    Eigen::Vector6d calculateDisplacementState( const EarthDeformationEpoch& epoch, const Eigen::Vector3d& position ) override;
    Eigen::Vector6d calculateDisplacementState( const EarthDeformationEpoch& epoch,
                                                const std::shared_ptr< ground_stations::GroundStationState > station ) override;
    //! Constructor creating pole tide deformation object from polar motion model.
    /*!
     *  Constructor creating pole tide deformation object from polar motion model.
     *  \param polarMotionCalculator Object capable of calculating pole position in ITRS at a given time instant.
     */
    PoleTideDeformation( const std::shared_ptr< earth_orientation::PolarMotionCalculator > polarMotionCalculator,
                         const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > contextClock = nullptr ):
        polarMotionCalculator_( polarMotionCalculator ), contextClock_( contextClock )
    {}

    ~PoleTideDeformation( ) {}

    //! Calculate displacement in the station's topocentric east, north, up frame, in metres.
    /*!
     *  Up follows the ellipsoid normal for an oblate station shape and the radial direction
     *  for a spherical or unspecified shape. The spherical IERS displacement is projected
     *  into this frame, so applying the inverse station topocentric rotation recovers
     *  calculateDisplacement for any supported shape.
     *  \param ephemerisTime Terrestrial time (TT), in seconds since J2000.
     *  \param nominalSiteState Nominal station position and topocentric frame.
     */
    Eigen::Vector3d calculateEnuDisplacement( const double ephemerisTime,
                                              const std::shared_ptr< ground_stations::GroundStationState > nominalSiteState );

    //! Calculate displacement in body-fixed Cartesian (ITRS) coordinates, in metres.
    /*!
     *  \param ephemerisTime Terrestrial time (TT), in seconds since J2000.
     *  \param stationState Nominal station position. The result uses the geocentric local
     *  basis and is independent of the shape used to represent this same Cartesian point.
     */
    Eigen::Vector3d calculateDisplacement( const double ephemerisTime,
                                           const std::shared_ptr< ground_stations::GroundStationState > stationState );

    Eigen::Vector3d calculateDisplacement( const double time, const Eigen::Vector3d& bodyFixedPosition )
    {
        return calculateDisplacement(
                time,
                std::make_shared< ground_stations::GroundStationState >( bodyFixedPosition, coordinate_conversions::cartesian_position ) );
    }

private:
    //! IERS Eq. (26) displacement in spherical east, north, radial-up coordinates, in metres.
    Eigen::Vector3d calculateGeocentricEnuDisplacement( const double ephemerisTime,
                                                        const std::shared_ptr< ground_stations::GroundStationState > nominalSiteState );

    //! Pointer to object to calculate the polar motion at a given time.
    /*!
     *  Pointer to object to calculate the polar motion at a given time.
     */
    std::shared_ptr< earth_orientation::PolarMotionCalculator > polarMotionCalculator_;
    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > contextClock_;
};

}  // namespace basic_astrodynamics

}  // namespace tudat
#endif  // TUDAT_POLETIDEDEFORMATION_H
