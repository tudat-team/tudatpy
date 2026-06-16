/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_FOURTHDEGREEFULLTWOBODYGRAVITATIONALTORQUEPARTIAL_H
#define TUDAT_FOURTHDEGREEFULLTWOBODYGRAVITATIONALTORQUEPARTIAL_H

#include "tudat/astro/gravitation/fourthDegreeFullTwoBodyGravitationalTorque.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityField.h"
#include "tudat/astro/orbit_determination/rotational_dynamics_partials/inertiaTensorPartial.h"
#include "tudat/astro/orbit_determination/rotational_dynamics_partials/torquePartial.h"

namespace tudat
{

namespace acceleration_partials
{

namespace detail
{

struct FourthDegreeTorqueAuxiliaryQuantities {
    double xCoordinate;
    double yCoordinate;
    double zCoordinate;
    double xCoordinateSquared;
    double yCoordinateSquared;
    double zCoordinateSquared;
    double xyTerm;
    double xzTerm;
    double yzTerm;
    double relativeDistanceSquared;
    double inverseRelativeDistanceSquared;
    double inverseRelativeDistanceToFourthPower;
    double relativeDistanceToFifthPower;
    double torquePrefactor;
    double aComponentOfBodyExertingTorque;
    double bComponentOfBodyExertingTorque;
    double cComponentOfBodyExertingTorque;
    double ixyComponentOfBodyExertingTorque;
    double ixzComponentOfBodyExertingTorque;
    double iyzComponentOfBodyExertingTorque;
    double traceOfInertiaTensorOfBodyExertingTorque;
    double contractedInertiaTensorOfBodyExertingTorque;
    double wPrimeQuantity;
    double fyzFunction;
    double fxzFunction;
    double fxyFunction;
    double gyzFunction;
    double gxzFunction;
    double gxyFunction;
};

Eigen::Matrix< double, 6, 1 > getIndependentInertiaTensorComponentsFromMatrix( const Eigen::Matrix3d& inertiaTensor );

FourthDegreeTorqueAuxiliaryQuantities computeFourthDegreeTorqueAuxiliaryQuantities(
        const Eigen::Vector3d& relativePositionOfBodyExertingTorqueInBodyFixedFrameOfBodyUndergoingTorque,
        const double massOfBodyExertingTorque,
        const Eigen::Matrix< double, 6, 1 >& independentInertiaTensorComponentsOfBodyExertingTorqueInFrameOfBodyUndergoingTorque );

double computePartialOfContractedInertiaTensorOfBodyExertingTorqueWrtCoordinate(
        const FourthDegreeTorqueAuxiliaryQuantities& auxiliaryQuantities,
        const int coordinateIndex );

Eigen::Matrix< double, 6, 1 > computePartialOfAuxiliaryFunctionsWrtPositionCoordinate(
        const FourthDegreeTorqueAuxiliaryQuantities& auxiliaryQuantities,
        const int coordinateIndex );

Eigen::Matrix< double, 6, 6 > computePartialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque(
        const FourthDegreeTorqueAuxiliaryQuantities& auxiliaryQuantities );

}  // namespace detail

//! Class to calculate partials of fourth-degree full two-body gravitational torque.
class FourthDegreeFullTwoBodyGravitationalTorquePartial : public TorquePartial
{
public:
    FourthDegreeFullTwoBodyGravitationalTorquePartial(
            const std::shared_ptr< gravitation::FourthDegreeFullTwoBodyGravitationalTorqueModel > torqueModel,
            const std::shared_ptr< gravitation::SphericalHarmonicsGravityField > gravityFieldOfBodyUndergoingTorque,
            const std::shared_ptr< gravitation::SphericalHarmonicsGravityField > gravityFieldOfBodyExertingTorque,
            const std::string& acceleratedBody,
            const std::string& acceleratingBody );

    ~FourthDegreeFullTwoBodyGravitationalTorquePartial( ) {}

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter ) override;

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter ) override;

    void wrtOrientationOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                          const bool addContribution = 1,
                                          const int startRow = 0,
                                          const int startColumn = 0 ) override;

    void wrtOrientationOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                           const bool addContribution = 1,
                                           const int startRow = 0,
                                           const int startColumn = 0 ) override;

    bool isStateDerivativeDependentOnIntegratedAdditionalStateTypes( const std::pair< std::string, std::string >& stateReferencePoint,
                                                                     const propagators::IntegratedStateType integratedStateType ) override;

    void wrtNonRotationalStateOfAdditionalBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                                const std::pair< std::string, std::string >& stateReferencePoint,
                                                const propagators::IntegratedStateType integratedStateType ) override;

    void update( const double currentTime = TUDAT_NAN ) override;

private:
    Eigen::Matrix< double, 6, 1 > getIndependentInertiaTensorComponentsFromMatrixDerivative(
            const Eigen::Matrix3d& inertiaTensorDerivative ) const;

    Eigen::Matrix3d getInertiaTensorPartialWrtNormalizedCosineCoefficient(
            const unsigned int order,
            const std::shared_ptr< gravitation::SphericalHarmonicsGravityField >& gravityField ) const;

    Eigen::Matrix3d getInertiaTensorPartialWrtNormalizedSineCoefficient(
            const unsigned int order,
            const std::shared_ptr< gravitation::SphericalHarmonicsGravityField >& gravityField ) const;

    Eigen::Matrix< double, 6, 1 > getIndependentInertiaTensorComponentsPartialWrtNormalizedCosineCoefficient(
            const unsigned int order,
            const std::shared_ptr< gravitation::SphericalHarmonicsGravityField >& gravityField ) const;

    Eigen::Matrix< double, 6, 1 > getIndependentInertiaTensorComponentsPartialWrtNormalizedSineCoefficient(
            const unsigned int order,
            const std::shared_ptr< gravitation::SphericalHarmonicsGravityField >& gravityField ) const;

    void wrtCosineSphericalHarmonicCoefficientsOfBodyUndergoingTorque( Eigen::MatrixXd& partialMatrix,
                                                                       const int c20Index,
                                                                       const int c21Index,
                                                                       const int c22Index );

    void wrtSineSphericalHarmonicCoefficientsOfBodyUndergoingTorque( Eigen::MatrixXd& partialMatrix,
                                                                     const int s21Index,
                                                                     const int s22Index );

    void wrtCosineSphericalHarmonicCoefficientsOfBodyExertingTorque( Eigen::MatrixXd& partialMatrix,
                                                                     const int c20Index,
                                                                     const int c21Index,
                                                                     const int c22Index );

    void wrtSineSphericalHarmonicCoefficientsOfBodyExertingTorque( Eigen::MatrixXd& partialMatrix, const int s21Index, const int s22Index );

    std::shared_ptr< gravitation::FourthDegreeFullTwoBodyGravitationalTorqueModel > torqueModel_;
    std::shared_ptr< gravitation::SphericalHarmonicsGravityField > gravityFieldOfBodyUndergoingTorque_;
    std::shared_ptr< gravitation::SphericalHarmonicsGravityField > gravityFieldOfBodyExertingTorque_;

    Eigen::Matrix< double, 3, 4 > currentPartialWrtQuaternionOfBodyUndergoingTorque_;
    Eigen::Matrix< double, 3, 4 > currentPartialWrtQuaternionOfBodyExertingTorque_;
    Eigen::Matrix3d currentPartialWrtPositionOfBodyUndergoingTorque_;
    Eigen::Matrix3d currentPartialWrtPositionOfBodyExertingTorque_;
    Eigen::Matrix< double, 3, 6 > currentPartialWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque_;
    Eigen::Matrix< double, 3, 6 > currentPartialWrtIndependentInertiaTensorComponentsOfBodyExertingTorqueInFrameOfBodyUndergoingTorque_;
    Eigen::Matrix3d currentRotationFromInertialToBodyFixedFrameOfBodyUndergoingTorque_;
    Eigen::Matrix3d currentRotationFromInertialToBodyFixedFrameOfBodyExertingTorque_;
    Eigen::Matrix3d currentRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorque_;
    Eigen::Vector3d currentRelativePositionOfBodyExertingTorqueInInertialFrame_;
    Eigen::Vector3d currentRelativePositionOfBodyExertingTorqueInBodyFixedFrameOfBodyUndergoingTorque_;
    Eigen::Matrix3d currentInertiaTensorOfBodyExertingTorque_;
};

}  // namespace acceleration_partials

}  // namespace tudat

#endif  // TUDAT_FOURTHDEGREEFULLTWOBODYGRAVITATIONALTORQUEPARTIAL_H
