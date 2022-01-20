#include <iostream>

#include "tudat/basics/utilities.h"
#include "tudat/math/basic/rotationRepresentations.h"

namespace tudat
{

namespace basic_mathematics
{

Eigen::Matrix< double, 4, 3 > calculateQuaternionWrtEulerAngle313Partial(
        const Eigen::Quaterniond& quaternion )
{
    double phiPlus = -std::atan2( quaternion.z( ), quaternion.w( ) );
    double cosPhiPlus = std::cos( phiPlus );
    double sinPhiPlus = std::sin( phiPlus );

    double phiMinus = -std::atan2( quaternion.y( ), -quaternion.x( ) );
    double cosPhiMinus = std::cos( phiMinus );
    double sinPhiMinus = std::sin( phiMinus );

    double cosineTheta = -quaternion.x( ) * quaternion.x( ) - quaternion.y( ) * quaternion.y( ) +
            quaternion.z( ) * quaternion.z( ) + quaternion.w( ) * quaternion.w( );
    if( cosineTheta > 1.0 )
    {
        cosineTheta = 1.0;
    }
    else if( cosineTheta < -1.0 )
    {
        cosineTheta = -1.0;
    }
    double thetaZero = std::acos( cosineTheta ) / 2.0;
    double cosThetaZero = std::cos( thetaZero );
    double sinThetaZero = std::sin( thetaZero );

    Eigen::MatrixXd scaledPartials = Eigen::MatrixXd::Zero( 4, 3 );

    scaledPartials << -sinThetaZero * cosPhiPlus, -cosThetaZero * sinPhiPlus, 0.0,
             -cosThetaZero * cosPhiMinus, 0.0, -sinThetaZero * sinPhiMinus,
              -cosThetaZero * sinPhiMinus, 0.0, sinThetaZero * cosPhiMinus,
              sinThetaZero * sinPhiPlus, -cosThetaZero * cosPhiPlus, 0.0;

    static Eigen::Matrix3d transformationMatrix =
            0.5 * ( Eigen::Matrix3d( ) << 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, -1.0 ).finished( );

    return scaledPartials * transformationMatrix;
}

//! Function to compute the partial derivative of 3-1-3 Euler angles w.r.t. entries of associated quaternion
Eigen::Matrix< double, 3, 4 > calculateEulerAngle313WrtQuaternionPartial(
        const Eigen::Quaterniond& quaternion )
{
    double termsSquared1 = ( quaternion.z( ) * quaternion.z( ) + quaternion.w( ) * quaternion.w( ) );
    double termsSquared2 = ( quaternion.x( ) * quaternion.x( ) + quaternion.y( ) * quaternion.y( ) );

    double recurrentTerm = std::sqrt( termsSquared2 / termsSquared1 );
    Eigen::Matrix< double, 3, 4 > partialMatrix;

    partialMatrix( 0, 0 ) = quaternion.z( ) / termsSquared1;
    partialMatrix( 0, 1 ) = quaternion.y( ) / termsSquared2;
    partialMatrix( 0, 2 ) = -quaternion.x( ) / termsSquared2;
    partialMatrix( 0, 3 ) = -quaternion.w( ) / termsSquared1;

    partialMatrix( 1, 0 ) = -2.0 * quaternion.w( ) * recurrentTerm;
    partialMatrix( 1, 1 ) = 2.0 * quaternion.x( ) / recurrentTerm;
    partialMatrix( 1, 2 ) = 2.0 * quaternion.y( ) / recurrentTerm;
    partialMatrix( 1, 3 ) = -2.0 * quaternion.z( ) * recurrentTerm;

    partialMatrix( 2, 0 ) = partialMatrix( 0, 0 );
    partialMatrix( 2, 1 ) = -partialMatrix( 0, 1 );
    partialMatrix( 2, 2 ) = -partialMatrix( 0, 2 );
    partialMatrix( 2, 3 ) = partialMatrix( 0, 3 );

    return partialMatrix;

}

//! Function to compute the partial derivative of 3-1-3 Euler angles w.r.t. entries of associated quaternion
Eigen::Matrix< double, 3, 4 > calculateEulerAngle313WrtQuaternionPartialFromEulerAngles(
        const Eigen::Vector3d& eulerAngles )
{
    return calculateEulerAngle313WrtQuaternionPartial(
                Eigen::Quaterniond(
                    Eigen::AngleAxisd( -eulerAngles( 0 ), Eigen::Vector3d::UnitZ( ) ) *
                    Eigen::AngleAxisd( -eulerAngles( 1 ), Eigen::Vector3d::UnitX( ) ) *
                    Eigen::AngleAxisd( -eulerAngles( 2 ), Eigen::Vector3d::UnitZ( ) ) ) );
}

//! Get quaternion from associated 3-1-3 Euler angles
Eigen::Quaterniond getQuaternionFrom313EulerAngles(
        const Eigen::Vector3d& eulerAngles )
{
    double cosineHalfTheta = std::cos( eulerAngles( 1 ) / 2.0 );
    double sineHalfTheta = std::sin( eulerAngles( 1 ) / 2.0 );

    return Eigen::Quaterniond( cosineHalfTheta * std::cos( ( eulerAngles( 0 ) + eulerAngles( 2 ) ) / 2.0 ),
                               -sineHalfTheta * std::cos( ( eulerAngles( 0 ) - eulerAngles( 2 ) ) / 2.0 ),
                               sineHalfTheta * std::sin( ( eulerAngles( 0 ) - eulerAngles( 2 ) ) / 2.0 ),
                               -cosineHalfTheta * std::sin( ( eulerAngles( 0 ) + eulerAngles( 2 ) ) / 2.0 ) );

}

//! Get classical 1-3-2 Euler angles set from rotation matrix
Eigen::Vector3d get132EulerAnglesFromRotationMatrix(
        const Eigen::Matrix3d& rotationMatrix )
{
    Eigen::Vector3d eulerAngles;
    eulerAngles( 0 ) = std::atan2( -rotationMatrix( 2, 1 ), rotationMatrix( 1, 1 ) );
    eulerAngles( 1 ) = std::asin( rotationMatrix( 0, 1 ) );
    eulerAngles( 2 ) = std::atan2( -rotationMatrix( 0, 2 ), rotationMatrix( 0, 0 ) );
    return eulerAngles;
}

//! Get classical 3-1-3 Euler angles set from quaternion
Eigen::Vector3d get313EulerAnglesFromQuaternion(
        const Eigen::Quaterniond& quaternion )
{

    double phiPlus = -std::atan2( quaternion.z( ), quaternion.w( ) );
    double phiMinus = -std::atan2( quaternion.y( ), -quaternion.x( ) );

    double psi = phiPlus - phiMinus;
    double phi = phiPlus + phiMinus;

    double cosineTheta = -quaternion.x( ) * quaternion.x( ) - quaternion.y( ) * quaternion.y( ) +
            quaternion.z( ) * quaternion.z( ) + quaternion.w( ) * quaternion.w( );

    double theta = std::acos( !( std::fabs( cosineTheta ) > 1.0 ) ? cosineTheta : utilities::sgn( cosineTheta ) * 1.0 );

    return ( Eigen::Vector3d( )<< psi, theta, phi ).finished( );
}

//! Get classical 3-1-3 Euler angles set from rotation matrix
Eigen::Vector3d get313EulerAnglesFromRotationMatrix(
        const Eigen::Matrix3d& rotationMatrix )
{
    double theta = std::acos( rotationMatrix( 2, 2 ) );
    double psi = std::atan2( rotationMatrix( 0, 2 ), rotationMatrix( 1, 2 ) );
    double phi = std::atan2( rotationMatrix( 2, 0 ), -rotationMatrix( 2, 1 ) );

    return ( Eigen::Vector3d( )<< psi, theta, phi ).finished( );
}


}

}
