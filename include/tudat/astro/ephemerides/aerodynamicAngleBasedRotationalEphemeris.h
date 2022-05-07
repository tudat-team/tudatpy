#ifndef TUDAT_AERODYNAMIC_ANGLE_BASED_ROTATIONAL_EPHEMERIS_H
#define TUDAT_AERODYNAMIC_ANGLE_BASED_ROTATIONAL_EPHEMERIS_H

#include "tudat/astro/ephemerides/rotationalEphemeris.h"

namespace tudat
{

namespace ephemerides
{

class AerodynamicAngleBasedRotationalEphemeris: public ephemerides::RotationalEphemeris
{
public:

    //! Constructor.
    /*!
     * Constructor, sets frames between which rotation is determined.
     * \param baseFrameOrientation Base frame identifier.
     * \param targetFrameOrientation Target frame identifier.
     */
    AerodynamicAngleBasedRotationalEphemeris(
            const std::function< Eigen::Vector3d( const double ) > aerodynamicAngleFunction,
            const std::string& baseFrameOrientation,
            const std::string& targetFrameOrientation, )
        : RotationalEphemeris( baseFrameOrientation, targetFrameOrientation )
    { }

    //! Virtual destructor.
    /*!
     * Virtual destructor.
     */
    virtual ~AerodynamicAngleBasedRotationalEphemeris( ) { }

    virtual Eigen::Quaterniond getRotationToBaseFrame(
            const double secondsSinceEpoch ) = 0;

    virtual Eigen::Quaterniond getRotationToTargetFrame(
            const double secondsSinceEpoch ) = 0;

    virtual Eigen::Matrix3d getDerivativeOfRotationToBaseFrame(
            const double secondsSinceEpoch ) = 0;

    virtual Eigen::Matrix3d getDerivativeOfRotationToTargetFrame(
            const double secondsSinceEpoch ) = 0;

protected:


    std::function< Eigen::Vector3d( const double ) > aerodynamicAngleFunction_;

};


} // namespace ephemerides

} // namespace tudat

#endif // TUDAT_AERODYNAMIC_ANGLE_BASED_ROTATIONAL_EPHEMERIS_H
