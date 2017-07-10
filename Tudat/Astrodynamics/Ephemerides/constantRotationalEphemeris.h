#ifndef CONSTANTROTATIONALEPHEMERIS_H
#define CONSTANTROTATIONALEPHEMERIS_H

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
    ConstantRotationalEphemeris( boost::function< Eigen::Quaterniond( ) > constantOrientationFunction,
                                 std::string originalFrameOrientation = "ECLIPJ2000",
                                 std::string targetFrameOrientation = "" ):
        RotationalEphemeris( originalFrameOrientation, targetFrameOrientation ),
        constantOrientationFunction_( constantOrientationFunction ) { }

    ConstantRotationalEphemeris( Eigen::Quaterniond constantOrientation,
                                 std::string originalFrameOrientation = "ECLIPJ2000",
                                 std::string targetFrameOrientation = "" ):
        RotationalEphemeris( originalFrameOrientation, targetFrameOrientation ),
        constantOrientationFunction_( boost::lambda::constant( constantOrientation ) ){ }

    Eigen::Quaterniond getRotationToBaseFrame( const double ephemerisTime )
    {
        return constantOrientationFunction_( );
    }

    Eigen::Matrix3d getDerivativeOfRotationFromFrame( const double ephemerisTime )
    {
        return Eigen::Matrix3d::Zero( );
    }

    Eigen::Quaterniond getConstantOrientation( )
    {
        return constantOrientationFunction_( );
    }

private:

    boost::function< Eigen::Quaterniond( ) > constantOrientationFunction_;

};

class WrapperRotationalEphemeris : public RotationalEphemeris
{
public:
    WrapperRotationalEphemeris( const boost::function< Eigen::Quaterniond( const double ) >& orientationFunction,
                                const std::string& originalFrameOrientation = "ECLIPJ2000",
                                const std::string& targetFrameOrientation = "" ):
        RotationalEphemeris( originalFrameOrientation, targetFrameOrientation ),
        orientationFunction_( orientationFunction ) { }

    Eigen::Quaterniond getRotationToBaseFrame( const double ephemerisTime )
    {
        return orientationFunction_( ephemerisTime );
    }

    Eigen::Matrix3d getDerivativeOfRotationFromFrame( const double ephemerisTime )
    {
        std::cerr<<"Error, derivative of rotation not implemented for wrapper rotation"<<std::endl;
        return Eigen::Matrix3d::Zero( );
    }

private:

    boost::function< Eigen::Quaterniond(  const double ) > orientationFunction_;

};


}

}

#endif // CONSTANTROTATIONALEPHEMERIS_H
