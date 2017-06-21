#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityField.h"
#include "Tudat/Astrodynamics/Propagators/rotationalMotionStateDerivative.h"


namespace tudat
{

namespace propagators
{

Eigen::Vector3d evaluateRotationalEquationsOfMotion(
        const Eigen::Matrix3d& inertiaTensor, const Eigen::Vector3d& totalTorque,
        const Eigen::Vector3d& rotationVector, const Eigen::Matrix3d& inertiaTimeDerivative )
{
    return inertiaTensor.inverse( ) * ( totalTorque - rotationVector.cross( inertiaTensor * rotationVector ) -
                                        inertiaTimeDerivative * rotationVector );
}

Eigen::Matrix4d getQuaterionToQuaternionRateMatrix( const Eigen::Vector3d& currentBodyFixedRotationRate )
{
    Eigen::Matrix4d conversionMatrix = Eigen::Matrix4d::Zero( );
    conversionMatrix( 1, 0 ) = currentBodyFixedRotationRate( 0 );
    conversionMatrix( 2, 0 ) = currentBodyFixedRotationRate( 1 );
    conversionMatrix( 3, 0 ) = currentBodyFixedRotationRate( 2 );

    conversionMatrix( 2, 1 ) = -currentBodyFixedRotationRate( 2 );
    conversionMatrix( 3, 1 ) = currentBodyFixedRotationRate( 1 );

    conversionMatrix( 3, 2 ) = -currentBodyFixedRotationRate( 0 );

    conversionMatrix( 0, 1 ) = -currentBodyFixedRotationRate( 0 );
    conversionMatrix( 0, 2 ) = -currentBodyFixedRotationRate( 1 );
    conversionMatrix( 0, 3 ) = -currentBodyFixedRotationRate( 2 );

    conversionMatrix( 1, 2 ) = currentBodyFixedRotationRate( 2 );
    conversionMatrix( 1, 3 ) = -currentBodyFixedRotationRate( 1 );

    conversionMatrix( 2, 3 ) = currentBodyFixedRotationRate( 0 );

   // std::cout<<"Matrix: "<<conversionMatrix<<std::endl;
    return 0.5 * ( conversionMatrix );
}


Eigen::Vector4d calculateQuaternionDerivative( const Eigen::Vector4d& currentQuaternion, const Eigen::Vector3d& currentBodyFixedRotationRate )
{
//    std::cout<<( getQuaterionToQuaternionRateMatrix( currentBodyFixedRotationRate ) * currentQuaternion ).transpose( )<<" "<<
//               currentBodyFixedRotationRate.norm( )<<" "<<currentQuaternion.transpose( )<<" "<<currentBodyFixedRotationRate.transpose( )<<std::endl;

//    std::cout<<"Rot: "<<currentBodyFixedRotationRate<<std::endl;
//    std::cout<<"Quat: "<<currentQuaternion<<std::endl;
//    std::cout<<"Der.: "<<getQuaterionToQuaternionRateMatrix( currentBodyFixedRotationRate ) * currentQuaternion<<std::endl;

//    sleep( 10000.0 );
    return getQuaterionToQuaternionRateMatrix( currentBodyFixedRotationRate ) * currentQuaternion;
}

}

}
