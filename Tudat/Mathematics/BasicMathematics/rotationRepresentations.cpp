#include <iostream>

#include "Tudat/Mathematics/BasicMathematics/rotationRepresentations.h"

namespace tudat
{

namespace basic_mathematics
{

Eigen::Matrix< double, 3, 4 > calculateEulerAngle313WrtQuaternionPartial(
        const Eigen::Quaterniond& quaternion )
{
    double termsSquared1 = ( quaternion.z( ) * quaternion.z( ) + quaternion.w( ) * quaternion.w( ) );
    double termsSquared2 = ( quaternion.x( ) * quaternion.x( ) + quaternion.y( ) * quaternion.y( ) );

    double recurrentTerm = std::sqrt( termsSquared2 / termsSquared1 );
    Eigen::Matrix< double, 3, 4 > partialMatrix;

    partialMatrix( 0, 0 ) = quaternion.z( ) / termsSquared1;
    partialMatrix( 0, 1 ) = -quaternion.y( ) / termsSquared2;
    partialMatrix( 0, 2 ) = quaternion.x( ) / termsSquared2;
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

Eigen::Matrix< double, 3, 4 > calculateEulerAngle313WrtQuaternionPartialFromEulerAngles(
        const Eigen::Vector3d& eulerAngles )
{
    return calculateEulerAngle313WrtQuaternionPartial(
                Eigen::Quaterniond(
                    Eigen::AngleAxisd( -eulerAngles( 2 ), Eigen::Vector3d::UnitZ( ) ) *
                    Eigen::AngleAxisd( -eulerAngles( 1 ), Eigen::Vector3d::UnitX( ) ) *
                    Eigen::AngleAxisd( -eulerAngles( 0 ), Eigen::Vector3d::UnitZ( ) ) ) );
}

Eigen::Quaterniond getQuaternionFrom313EulerAngles(
        const Eigen::Vector3d& eulerAngles )
{
    double cosineHalfTheta = std::cos( eulerAngles( 1 ) / 2.0 );
    double sineHalfTheta = std::sin( eulerAngles( 1 ) / 2.0 );

    return Eigen::Quaterniond( cosineHalfTheta * std::cos( ( eulerAngles( 0 ) + eulerAngles( 2 ) ) / 2.0 ),
                               -sineHalfTheta * std::cos( ( eulerAngles( 0 ) - eulerAngles( 2 ) ) / 2.0 ),
                               -sineHalfTheta * std::sin( ( eulerAngles( 0 ) - eulerAngles( 2 ) ) / 2.0 ),
                               -cosineHalfTheta * std::sin( ( eulerAngles( 0 ) + eulerAngles( 2 ) ) / 2.0 ) );

}

Eigen::Vector3d get313EulerAnglesFromQuaternion(
        const Eigen::Quaterniond& quaternion )
{
    double theta = 2.0 * atan2( std::sqrt( quaternion.x( ) * quaternion.x( ) + quaternion.y( ) * quaternion.y( ) ),
                                std::sqrt( quaternion.z( ) * quaternion.z( ) + quaternion.w( ) * quaternion.w( ) ) );
    double phiPlus = atan2( -quaternion.z( ), quaternion.w( ) );
    double phiMinus = atan2( -quaternion.y( ), -quaternion.x( ) );
    Eigen::Vector3d eulerAngles = ( Eigen::Vector3d( )<< phiPlus + phiMinus, theta, phiPlus - phiMinus ).finished( );

    return eulerAngles;

}

}

}
