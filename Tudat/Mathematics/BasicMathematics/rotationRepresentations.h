#ifndef TUDAT_ROTATIONREPRESENTATIONS_H
#define TUDAT_ROTATIONREPRESENTATIONS_H

#include <Eigen/Core>
#include <Eigen/Dense>

namespace tudat
{

namespace basic_mathematics
{

//Eigen::Vector4d getManualQuaternionEntriesFromRotationMatrix(
//        const Eigen::Matrix3d rotationMatrix );

//Eigen::Matrix3d getManualRotationMatrixFromQuaternion(
//        const Eigen::Vector4d rotationMatrix );

//class RotationMatrixPartialCalculator
//{
//public:

//    RotationMatrixPartialCalculator( )
//    {
//        matrixA_.setZero( );
//        partial13Multiplier_.setZero( );
//        partial23Multiplier_.setZero( );
//        partial31Multiplier_.setZero( );

//    }

//    void calculateRotationMatrixPartialDerivatives(
//            const Eigen::Matrix3d& currentRotationMatrix,
//            const std::vector< Eigen::Vector3d > independentPartials,
//            Eigen::Matrix< double, 6, Eigen::Dynamic >& dependentPartialsMatrix );

//private:

//    void setMatrixA( const Eigen::Matrix3d& currentRotationMatrix );

//    void setIndependentPartialMultipliers( const Eigen::Matrix3d& currentRotationMatrix );

//    Eigen::Matrix< double, 6, 6 > matrixA_;

//    Eigen::Matrix< double, 6, 1 > partial13Multiplier_;

//    Eigen::Matrix< double, 6, 1 > partial23Multiplier_;

//    Eigen::Matrix< double, 6, 1 > partial31Multiplier_;

//    Eigen::Matrix< double, 6, Eigen::Dynamic > rightHandSide_;

//};

//double calculateDependentQuaternionPartial(
//        const Eigen::Quaterniond& currentQuaternion,
//        const Eigen::Vector3d independentEntryPartials );

//void calculateFullQuaternionPartials(
//        const Eigen::Quaterniond& currentQuaternion,
//        const Eigen::Vector3d independentEntryPartials,
//        Eigen::Vector4d& fullPartial );

//void calculatePartialOfQuaternionWrtRotationMatrix(
//        const Eigen::Vector4d rotationMatrix,
//        std::vector< Eigen::Matrix3d >& partials );

//std::vector< Eigen::Vector3d > getPartialsOfIndependentRotationMatrixEntriesWrtQuaternion(
//            Eigen::Vector4d quaternionEntries );

//void calculatePartialOfRotationMatrixWrtQuaternion(
//        const Eigen::Matri/x3d rotationMatrix,
//        std::vector< Eigen::Matrix3d >& partials );

Eigen::Matrix< double, 3, 4 > calculateEulerAngle313WrtQuaternionPartial(
        const Eigen::Quaterniond& quaternion );

Eigen::Matrix< double, 3, 4 > calculateEulerAngle313WrtQuaternionPartialFromEulerAngles(
        const Eigen::Vector3d& eulerAngles );

Eigen::Quaterniond getQuaternionFrom313EulerAngles(
        const Eigen::Vector3d& eulerAngles );

Eigen::Vector3d get313EulerAnglesFromQuaternion(
        const Eigen::Quaterniond& quaternion );
}

}

#endif // TUDAT_ROTATIONREPRESENTATIONS_H
