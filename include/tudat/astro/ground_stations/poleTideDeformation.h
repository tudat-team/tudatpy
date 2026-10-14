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
    PoleTideDeformation( const std::shared_ptr< earth_orientation::PolarMotionCalculator > polarMotionCalculator,
                         const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > contextClock = nullptr ):
        BodyDeformationModel( contextClock ), polarMotionCalculator_( polarMotionCalculator )
    {}

    ~PoleTideDeformation( ) {}

protected:
    Eigen::Vector3d evaluateDisplacement( const Time& tdb, const std::shared_ptr< ground_stations::GroundStationState > station ) override;
    Eigen::Vector6d evaluateDisplacementState( const Time& tdb,
                                               const std::shared_ptr< ground_stations::GroundStationState > station ) override;

private:
    Eigen::Vector3d displacementAtTt( const Time& tt, const std::shared_ptr< ground_stations::GroundStationState > station );

    //! Pointer to object to calculate the polar motion at a given time.
    /*!
     *  Pointer to object to calculate the polar motion at a given time.
     */
    std::shared_ptr< earth_orientation::PolarMotionCalculator > polarMotionCalculator_;
};

}  // namespace basic_astrodynamics

}  // namespace tudat
#endif  // TUDAT_POLETIDEDEFORMATION_H
