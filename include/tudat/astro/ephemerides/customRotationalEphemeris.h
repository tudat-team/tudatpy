#ifndef TUDAT_CUSTOM_ROTATIONAL_EPHEMERIS_H
#define TUDAT_CUSTOM_ROTATIONAL_EPHEMERIS_H

#include "tudat/astro/ephemerides/rotationalEphemeris.h"

namespace tudat
{

namespace ephemerides
{

class CustomRotationalEphemeris: public ephemerides::RotationalEphemeris
{
public:

    //! Constructor.
    /*!
     * Constructor, sets frames between which rotation is determined.
     * \param baseFrameOrientation Base frame identifier.
     * \param targetFrameOrientation Target frame identifier.
     */
    CustomRotationalEphemeris(
            const std::function< Eigen::Quaterniond( const double ) > targetToBaseFrameOrientationFunction,
            const std::string& baseFrameOrientation,
            const std::string& targetFrameOrientation,
            const double numericalDerivativeTimeStep = 10.0 )
        : RotationalEphemeris( baseFrameOrientation, targetFrameOrientation ),
          targetToBaseFrameOrientationFunction_( targetToBaseFrameOrientationFunction ),
          numericalDerivativeTimeStep_( numericalDerivativeTimeStep )
    { }

    //! Virtual destructor.
    /*!
     * Virtual destructor.
     */
    virtual ~CustomRotationalEphemeris( ) { }

    virtual Eigen::Quaterniond getRotationToBaseFrame(
            const double secondsSinceEpoch )
    {
        return targetToBaseFrameOrientationFunction_( secondsSinceEpoch );
    }

    virtual Eigen::Quaterniond getRotationToTargetFrame(
            const double secondsSinceEpoch )
    {
        return getRotationToBaseFrame( secondsSinceEpoch ).inverse( );
    }

    virtual Eigen::Matrix3d getDerivativeOfRotationToBaseFrame(
            const double secondsSinceEpoch )
    {
        return ( targetToBaseFrameOrientationFunction_( secondsSinceEpoch + numericalDerivativeTimeStep_ ).toRotationMatrix( ) -
                 targetToBaseFrameOrientationFunction_( secondsSinceEpoch - numericalDerivativeTimeStep_ ).toRotationMatrix( ) ) /
                ( 2.0 * numericalDerivativeTimeStep_ );
    }

    virtual Eigen::Matrix3d getDerivativeOfRotationToTargetFrame(
            const double secondsSinceEpoch )
    {
        return getDerivativeOfRotationToBaseFrame( secondsSinceEpoch ).transpose( );
    }

protected:

    std::function< Eigen::Quaterniond( const double ) > targetToBaseFrameOrientationFunction_;

    double numericalDerivativeTimeStep_;
};


} // namespace ephemerides

} // namespace tudat

#endif // TUDAT_CUSTOM_ROTATIONAL_EPHEMERIS_H
