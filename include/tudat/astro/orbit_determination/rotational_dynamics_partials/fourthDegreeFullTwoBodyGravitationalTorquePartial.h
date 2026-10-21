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
#include "tudat/astro/orbit_determination/observation_partials/rotationMatrixPartial.h"
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
            const std::string& acceleratingBody,
            const observation_partials::RotationMatrixPartialNamedList& rotationMatrixPartialsOfBodyUndergoingTorque =
                    observation_partials::RotationMatrixPartialNamedList( ),
            const observation_partials::RotationMatrixPartialNamedList& rotationMatrixPartialsOfBodyExertingTorque =
                    observation_partials::RotationMatrixPartialNamedList( ) );

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

    void resetCurrentTimeOfMemberObjects( ) override
    {
        torqueModel_->resetCurrentTime( );
    }

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

    void wrtCosineSphericalHarmonicCoefficientsOfBodyUndergoingTorque( Eigen::MatrixXd& partialMatrix,
                                                                       const std::vector< std::pair< int, int > >& blockIndices );

    void wrtSineSphericalHarmonicCoefficientsOfBodyUndergoingTorque( Eigen::MatrixXd& partialMatrix,
                                                                     const int s21Index,
                                                                     const int s22Index );

    void wrtSineSphericalHarmonicCoefficientsOfBodyUndergoingTorque( Eigen::MatrixXd& partialMatrix,
                                                                     const std::vector< std::pair< int, int > >& blockIndices );

    void wrtCosineSphericalHarmonicCoefficientsOfBodyExertingTorque( Eigen::MatrixXd& partialMatrix,
                                                                     const int c20Index,
                                                                     const int c21Index,
                                                                     const int c22Index );

    void wrtCosineSphericalHarmonicCoefficientsOfBodyExertingTorque( Eigen::MatrixXd& partialMatrix,
                                                                     const std::vector< std::pair< int, int > >& blockIndices );

    void wrtSineSphericalHarmonicCoefficientsOfBodyExertingTorque( Eigen::MatrixXd& partialMatrix, const int s21Index, const int s22Index );

    void wrtSineSphericalHarmonicCoefficientsOfBodyExertingTorque( Eigen::MatrixXd& partialMatrix,
                                                                   const std::vector< std::pair< int, int > >& blockIndices );

    void wrtMassOfBodyExertingTorque( Eigen::MatrixXd& partialMatrix );

    void wrtRotationModelParameter( Eigen::MatrixXd& partialMatrix,
                                    const bool wrtBodyUndergoingTorque,
                                    const estimatable_parameters::EstimatebleParametersEnum parameterType,
                                    const std::string& secondaryIdentifier );

    void wrtPolynomialGravityFieldVariations( const bool wrtBodyUndergoingTorque,
                                              const std::vector< std::pair< int, int > >& cosineBlockIndices,
                                              const std::vector< std::pair< int, int > >& sineBlockIndices,
                                              const std::vector< std::vector< std::pair< int, int > > > powersPerCosineBlockIndex,
                                              const std::vector< std::vector< std::pair< int, int > > > powersPerSineBlockIndex,
                                              const double referenceEpoch,
                                              Eigen::MatrixXd& partialMatrix );

    void wrtPeriodicGravityFieldVariations( const bool wrtBodyUndergoingTorque,
                                            const std::vector< std::pair< int, int > >& cosineBlockIndices,
                                            const std::vector< std::pair< int, int > >& sineBlockIndices,
                                            const std::vector< std::vector< std::pair< int, int > > > periodsPerCosineBlockIndex,
                                            const std::vector< std::vector< std::pair< int, int > > > periodsPerSineBlockIndex,
                                            const std::vector< double >& frequencies,
                                            const double referenceEpoch,
                                            Eigen::MatrixXd& partialMatrix );

    std::shared_ptr< gravitation::FourthDegreeFullTwoBodyGravitationalTorqueModel > torqueModel_;
    std::shared_ptr< gravitation::SphericalHarmonicsGravityField > gravityFieldOfBodyUndergoingTorque_;
    std::shared_ptr< gravitation::SphericalHarmonicsGravityField > gravityFieldOfBodyExertingTorque_;

    Eigen::Matrix< double, 3, 4 > currentPartialWrtQuaternionOfBodyUndergoingTorque_;
    Eigen::Matrix< double, 3, 4 > currentPartialWrtQuaternionOfBodyExertingTorque_;
    Eigen::Matrix3d currentPartialWrtPositionOfBodyUndergoingTorque_;
    Eigen::Matrix3d currentPartialWrtPositionOfBodyExertingTorque_;
    Eigen::Matrix< double, 3, 6 > currentPartialWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque_;
    Eigen::Matrix< double, 3, 6 > currentPartialWrtIndependentInertiaTensorComponentsOfBodyExertingTorqueInFrameOfBodyUndergoingTorque_;
    Eigen::Vector3d currentPartialWrtMassOfBodyExertingTorque_;
    Eigen::Matrix3d currentRotationFromInertialToBodyFixedFrameOfBodyUndergoingTorque_;
    Eigen::Matrix3d currentRotationFromInertialToBodyFixedFrameOfBodyExertingTorque_;
    Eigen::Matrix3d currentRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorque_;
    Eigen::Vector3d currentRelativePositionOfBodyExertingTorqueInInertialFrame_;
    Eigen::Vector3d currentRelativePositionOfBodyExertingTorqueInBodyFixedFrameOfBodyUndergoingTorque_;
    Eigen::Matrix3d currentInertiaTensorOfBodyExertingTorque_;

    observation_partials::RotationMatrixPartialNamedList rotationMatrixPartialsOfBodyUndergoingTorque_;
    observation_partials::RotationMatrixPartialNamedList rotationMatrixPartialsOfBodyExertingTorque_;
};

}  // namespace acceleration_partials

}  // namespace tudat

#endif  // TUDAT_FOURTHDEGREEFULLTWOBODYGRAVITATIONALTORQUEPARTIAL_H
