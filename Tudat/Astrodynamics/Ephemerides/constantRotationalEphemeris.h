#ifndef TUDAT_CONSTANTROTATIONALEPHEMERIS_H
#define TUDAT_CONSTANTROTATIONALEPHEMERIS_H

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"
#include "Tudat/Basics/utilities.h"

namespace tudat
{

namespace ephemerides
{

class ConstantRotationalEphemeris : public RotationalEphemeris
{
public:

    ConstantRotationalEphemeris( Eigen::Quaterniond constantOrientation,
                                 std::string originalFrameOrientation = "ECLIPJ2000",
                                 std::string targetFrameOrientation = "" ):
        RotationalEphemeris( originalFrameOrientation, targetFrameOrientation ),
        constantOrientation_( constantOrientation ){ }

    Eigen::Quaterniond getRotationToBaseFrame( const double ephemerisTime )
    {
        return constantOrientation_;
    }

    Eigen::Quaterniond getRotationToTargetFrame( const double ephemerisTime )
    {
        return constantOrientation_.inverse( );
    }

    Eigen::Matrix3d getDerivativeOfRotationToBaseFrame( const double ephemerisTime )
    {
        return Eigen::Matrix3d::Zero( );
    }

    Eigen::Matrix3d getDerivativeOfRotationToTargetFrame( const double ephemerisTime )
    {
        return Eigen::Matrix3d::Zero( );
    }

    Eigen::Quaterniond getConstantOrientation( )
    {
        return constantOrientation_;
    }

private:

    Eigen::Quaterniond constantOrientation_;

};

}

}

#endif // TUDAT_CONSTANTROTATIONALEPHEMERIS_H
