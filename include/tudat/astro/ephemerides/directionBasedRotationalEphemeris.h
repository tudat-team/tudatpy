#ifndef TUDAT_DIRECTION_BASED_ROTATIONAL_EPHEMERIS_H
#define TUDAT_DIRECTION_BASED_ROTATIONAL_EPHEMERIS_H

#include "tudat/astro/ephemerides/rotationalEphemeris.h"

namespace tudat
{

namespace ephemerides
{

class DirectionBasedRotationalEphemeris: public ephemerides::RotationalEphemeris
{
public:

    //! Constructor.
    /*!
     * Constructor, sets frames between which rotation is determined.
     * \param baseFrameOrientation Base frame identifier.
     * \param targetFrameOrientation Target frame identifier.
     */
    DirectionBasedRotationalEphemeris(
            const std::function< Eigen::Vector3d( const double ) > inertialBodyAxisDirectionFunction,
            const std::string& baseFrameOrientation,
            const std::string& targetFrameOrientation, )
        : RotationalEphemeris( baseFrameOrientation, targetFrameOrientation )
    { }

    //! Virtual destructor.
    /*!
     * Virtual destructor.
     */
    virtual ~DirectionBasedRotationalEphemeris( ) { }

    virtual Eigen::Quaterniond getRotationToBaseFrame(
            const double secondsSinceEpoch ) = 0;

    virtual Eigen::Quaterniond getRotationToTargetFrame(
            const double secondsSinceEpoch ) = 0;

    virtual Eigen::Matrix3d getDerivativeOfRotationToBaseFrame(
            const double secondsSinceEpoch ) = 0;

    virtual Eigen::Matrix3d getDerivativeOfRotationToTargetFrame(
            const double secondsSinceEpoch ) = 0;

protected:

    //! Function returning thrust-direction (represented in the relevant propagation frame) as a function of time.
    std::function< Eigen::Vector3d( const double ) > forceDirectionFunction_;
};


} // namespace ephemerides

} // namespace tudat

#endif // TUDAT_DIRECTION_BASED_ROTATIONAL_EPHEMERIS_H
