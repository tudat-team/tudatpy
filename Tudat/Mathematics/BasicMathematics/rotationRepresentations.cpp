#include <iostream>

#include "Tudat/Mathematics/BasicMathematics/rotationRepresentations.h"

namespace tudat
{

namespace basic_mathematics
{

//Eigen::Vector4d getManualQuaternionEntriesFromRotationMatrix(
//        const Eigen::Matrix3d rotationMatrix )
//{
//    Eigen::Vector4d vector;
//    vector( 0 ) = 0.5 * std::sqrt( 1.0 + rotationMatrix( 0, 0 ) + rotationMatrix( 1, 1 ) + rotationMatrix( 2, 2 ) );
//    vector( 1 ) = ( rotationMatrix( 2, 1 ) - rotationMatrix( 1, 2 ) ) / ( 4.0 * vector( 0 ) );
//    vector( 2 ) = ( rotationMatrix( 0, 2 ) - rotationMatrix( 2, 0 ) ) / ( 4.0 * vector( 0 ) );
//    vector( 3 ) = ( rotationMatrix( 1, 0 ) - rotationMatrix( 0, 1 ) ) / ( 4.0 * vector( 0 ) );

//    return vector;
//}

//Eigen::Matrix3d getManualRotationMatrixFromQuaternion(
//        const Eigen::Vector4d quaternionVector )
//{
//    Eigen::Quaterniond quaternion =
//            Eigen::Quaterniond( quaternionVector( 0 ), quaternionVector( 1 ), quaternionVector( 2 ), quaternionVector( 3 ) );
//    return quaternion.toRotationMatrix( );
//}


//void RotationMatrixPartialCalculator::calculateRotationMatrixPartialDerivatives(
//        const Eigen::Matrix3d& currentRotationMatrix,
//        const std::vector< Eigen::Vector3d > independentPartials,
//        Eigen::Matrix< double, 6, Eigen::Dynamic >& dependentPartialsMatrix )
//{
//    setMatrixA( currentRotationMatrix );
//    setIndependentPartialMultipliers( currentRotationMatrix );

//    rightHandSide_.setZero( 6, independentPartials.size( ) );
//    for( unsigned int i = 0; i < independentPartials.size( ); i++ )
//    {
//        rightHandSide_.block( 0, i, 6, 1 ) = independentPartials.at( i )( 0 ) * partial13Multiplier_ +
//                independentPartials.at( i )( 1 ) * partial23Multiplier_ +
//                independentPartials.at( i )( 2 ) * partial31Multiplier_;
//    }

//    std::cout<<"Mat: "<<matrixA_<<std::endl;

//    dependentPartialsMatrix = matrixA_.partialPivLu( ).solve( rightHandSide_ );
//}

//void RotationMatrixPartialCalculator::setMatrixA( const Eigen::Matrix3d& currentRotationMatrix )
//{
//    matrixA_( 0, 0 ) = currentRotationMatrix( 0, 0 );
//    matrixA_( 0, 1 ) = currentRotationMatrix( 1, 0 );

//    matrixA_( 1, 2 ) = currentRotationMatrix( 0, 1 );
//    matrixA_( 1, 3 ) = currentRotationMatrix( 1, 1 );

//    matrixA_( 2, 4 ) = currentRotationMatrix( 1, 2 );
//    matrixA_( 2, 5 ) = currentRotationMatrix( 2, 2 );

//    matrixA_( 3, 0 ) = currentRotationMatrix( 0, 1 );
//    matrixA_( 3, 1 ) = currentRotationMatrix( 1, 1 );
//    matrixA_( 3, 2 ) = currentRotationMatrix( 0, 0 );
//    matrixA_( 3, 3 ) = currentRotationMatrix( 1, 0 );

//    matrixA_( 4, 0 ) = currentRotationMatrix( 0, 2 );
//    matrixA_( 4, 1 ) = currentRotationMatrix( 1, 2 );
//    matrixA_( 4, 4 ) = currentRotationMatrix( 1, 0 );
//    matrixA_( 4, 5 ) = currentRotationMatrix( 2, 0 );

//    matrixA_( 5, 2 ) = currentRotationMatrix( 0, 2 );
//    matrixA_( 5, 3 ) = currentRotationMatrix( 1, 2 );
//    matrixA_( 5, 4 ) = currentRotationMatrix( 1, 1 );
//    matrixA_( 5, 5 ) = currentRotationMatrix( 2, 1 );
//}

//void RotationMatrixPartialCalculator::setIndependentPartialMultipliers( const Eigen::Matrix3d& currentRotationMatrix )
//{
//    partial13Multiplier_( 0 ) = -currentRotationMatrix( 2, 0 );
//    partial13Multiplier_( 3 ) = -currentRotationMatrix( 2, 1 );
//    partial13Multiplier_( 4 ) = -currentRotationMatrix( 2, 2 );

//    partial23Multiplier_( 1 ) = -currentRotationMatrix( 2, 1 );
//    partial23Multiplier_( 3 ) = -currentRotationMatrix( 2, 0 );
//    partial23Multiplier_( 5 ) = -currentRotationMatrix( 2, 2 );

//    partial31Multiplier_( 2 ) = -currentRotationMatrix( 0, 2 );
//    partial31Multiplier_( 4 ) = -currentRotationMatrix( 0, 0 );
//    partial31Multiplier_( 5 ) = -currentRotationMatrix( 0, 1 );

//}

//double calculateDependentQuaternionPartial(
//        const Eigen::Quaterniond& currentQuaternion,
//        const Eigen::Vector3d independentEntryPartials )
//{
//    return -( currentQuaternion.x( ) * independentEntryPartials.x( ) +
//              currentQuaternion.y( ) * independentEntryPartials.y( ) +
//              currentQuaternion.z( ) * independentEntryPartials.z( ) ) / currentQuaternion.w( );
//}

//void calculateFullQuaternionPartials(
//        const Eigen::Quaterniond& currentQuaternion,
//        const Eigen::Vector3d independentEntryPartials,
//        Eigen::Vector4d& fullPartial )
//{
//    fullPartial = ( Eigen::Vector4d( )<<calculateDependentQuaternionPartial(
//                        currentQuaternion, independentEntryPartials ), independentEntryPartials ).finished( );
//}

//void calculatePartialOfQuaternionWrtRotationMatrix(
//        const Eigen::Vector4d quaternion,
//        std::vector< Eigen::Matrix3d >& partials )
//{
//    partials.resize( 4 );
//    partials[ 0 ] = Eigen::Matrix3d::Identity( ) / ( 8.0 * quaternion( 0 ) );
//    for( unsigned int i = 1; i < 4; i++ )
//    {
//        partials[ i ] = -0.5 * Eigen::Matrix3d::Identity( ) * quaternion( i ) / quaternion( 0 );
//    }

//    partials[ 1 ]( 1, 2 ) += -1.0;
//    partials[ 1 ]( 2, 1 ) += 1.0;

//    partials[ 2 ]( 2, 0 ) += -1.0;
//    partials[ 2 ]( 0, 2 ) += 1.0;

//    partials[ 3 ]( 0, 1 ) += -1.0;
//    partials[ 3 ]( 1, 0 ) += 1.0;

//    for( unsigned int i = 1; i < 4; i++ )
//    {
//        partials[ i ] *= 1.0 / ( 4.0 * quaternion( 0 ) );
//    }
//}

//std::vector< Eigen::Vector3d > getPartialsOfIndependentRotationMatrixEntriesWrtQuaternion(
//            Eigen::Vector4d quaternionEntries )
//{
//     std::vector< Eigen::Vector3d > independentPartials;
//     independentPartials.push_back(
//                 2.0 * ( Eigen::Vector3d( )<< quaternionEntries( 2 ), -quaternionEntries( 1 ), -quaternionEntries( 2 ) ).finished( ) );
//     independentPartials.push_back(
//                 2.0 * ( Eigen::Vector3d( )<< quaternionEntries( 3 ), -quaternionEntries( 0 ), quaternionEntries( 3 ) ).finished( ) );
//     independentPartials.push_back(
//                 2.0 * ( Eigen::Vector3d( )<< quaternionEntries( 0 ), quaternionEntries( 3 ), -quaternionEntries( 0 ) ).finished( ) );
//     independentPartials.push_back(
//                 2.0 * ( Eigen::Vector3d( )<< quaternionEntries( 1 ), quaternionEntries( 2 ), quaternionEntries( 1 ) ).finished( ) );

//     std::cout<<"Indep partials"<<std::endl<<
//                independentPartials[ 0 ].transpose( )<<std::endl<<
//     independentPartials[ 1 ].transpose( )<<std::endl<<
//     independentPartials[ 2 ].transpose( )<<std::endl<<
//     independentPartials[ 3 ].transpose( )<<std::endl;

//     return independentPartials;
//}

//void calculatePartialOfRotationMatrixWrtQuaternion(
//        const Eigen::Matrix3d rotationMatrix,
//        std::vector< Eigen::Matrix3d >& partials )
//{
//    Eigen::Vector4d quaternionEntries = getManualQuaternionEntriesFromRotationMatrix(
//            rotationMatrix );
//    std::vector< Eigen::Vector3d > independentPartials = getPartialsOfIndependentRotationMatrixEntriesWrtQuaternion(
//                quaternionEntries );

//    RotationMatrixPartialCalculator rotationMatrixPartialCalculator;

//    Eigen::Matrix< double, 6, Eigen::Dynamic > dependentPartialsMatrix;

//    rotationMatrixPartialCalculator.calculateRotationMatrixPartialDerivatives(
//            rotationMatrix, independentPartials, dependentPartialsMatrix );

//    partials.resize( 4 );
//    for( int i = 0; i < 4; i++ )
//    {
//        partials[ i ]( 0, 0 ) = dependentPartialsMatrix( 0, i );
//        partials[ i ]( 0, 1 ) = dependentPartialsMatrix( 1, i );
//        partials[ i ]( 0, 2 ) = independentPartials[ i ]( 0 );
//        partials[ i ]( 1, 0 ) = dependentPartialsMatrix( 2, i );
//        partials[ i ]( 1, 1 ) = dependentPartialsMatrix( 3, i );
//        partials[ i ]( 1, 2 ) = independentPartials[ i ]( 1 );
//        partials[ i ]( 2, 0 ) = independentPartials[ i ]( 2 );
//        partials[ i ]( 2, 1 ) = dependentPartialsMatrix( 4, i );
//        partials[ i ]( 2, 2 ) = dependentPartialsMatrix( 5, i );
//        std::cout<<"Setting partials "<<i<<" "
//                <<independentPartials[ i ]( 0 )<<" "<<independentPartials[ i ]( 1 )<<" "<<independentPartials[ i ]( 2 )<<std::endl;
//    }
//}
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
