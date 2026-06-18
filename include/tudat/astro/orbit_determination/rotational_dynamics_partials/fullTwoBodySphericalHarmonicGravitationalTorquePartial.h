/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_FULLTWOBODYSPHERICALHARMONICGRAVITATIONALTORQUEPARTIAL_H
#define TUDAT_FULLTWOBODYSPHERICALHARMONICGRAVITATIONALTORQUEPARTIAL_H

#include "tudat/astro/gravitation/fullTwoBodySphericalHarmonicTorque.h"
#include "tudat/astro/orbit_determination/acceleration_partials/fullTwoBodySphericalHarmonicGravityPartial.h"
#include "tudat/astro/orbit_determination/rotational_dynamics_partials/torquePartial.h"

namespace tudat
{

namespace acceleration_partials
{

//! Class to calculate analytical partials of full two-body spherical-harmonic torque.
/*!
 * Provides derivatives of the full two-body torque model based on Dirkx et al. (2019):
 * body-2 spin torque from Dirkx et al. (2019), Eq. (60), with coefficient summation from Dirkx et al. (2019), Eq. (67),
 * and total/body torque relation from Dirkx et al. (2019), Eqs. (68)-(69), using effective coefficients from
 * Dirkx et al. (2019), Eqs. (47)-(49).
 */
class FullTwoBodySphericalHarmonicGravitationalTorquePartial : public TorquePartial
{
public:
    //! Constructor.
    FullTwoBodySphericalHarmonicGravitationalTorquePartial(
            const std::shared_ptr< gravitation::FullTwoBodySphericalHarmonicTorque > torqueModel,
            const std::shared_ptr< FullTwoBodySphericalHarmonicsGravityPartial > accelerationPartial,
            const std::string& acceleratedBody,
            const std::string& acceleratingBody,
            const observation_partials::RotationMatrixPartialNamedList& rotationMatrixPartialsOfBodyUndergoingTorque =
                    observation_partials::RotationMatrixPartialNamedList( ),
            const observation_partials::RotationMatrixPartialNamedList& rotationMatrixPartialsOfBodyExertingTorque =
                    observation_partials::RotationMatrixPartialNamedList( ) );

    ~FullTwoBodySphericalHarmonicGravitationalTorquePartial( ) {}

    //! Retrieve scalar-parameter partial function (none implemented).
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter ) override;

    //! Retrieve vector-parameter partial function (spherical-harmonic coefficient blocks).
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter ) override;

    //! Insert partial w.r.t. quaternion of body undergoing torque.
    void wrtOrientationOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                          const bool addContribution = 1,
                                          const int startRow = 0,
                                          const int startColumn = 0 ) override;

    //! Insert partial w.r.t. quaternion of body exerting torque.
    void wrtOrientationOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                           const bool addContribution = 1,
                                           const int startRow = 0,
                                           const int startColumn = 0 ) override;

    //! Return whether the torque derivative depends on a given additional integrated state type.
    bool isStateDerivativeDependentOnIntegratedAdditionalStateTypes( const std::pair< std::string, std::string >& stateReferencePoint,
                                                                     const propagators::IntegratedStateType integratedStateType ) override;

    //! Insert partial w.r.t. non-rotational state of an additional body.
    void wrtNonRotationalStateOfAdditionalBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                                const std::pair< std::string, std::string >& stateReferencePoint,
                                                const propagators::IntegratedStateType integratedStateType ) override;

    //! Update all cached torque partial terms to current model state.
    void update( const double currentTime = TUDAT_NAN ) override;

private:
    //! Partial w.r.t. cosine SH coefficient block of body undergoing torque.
    void wrtCosineSphericalHarmonicCoefficientsOfBodyUndergoingTorque( Eigen::MatrixXd& partialMatrix,
                                                                       const std::vector< std::pair< int, int > >& blockIndices );

    //! Partial w.r.t. sine SH coefficient block of body undergoing torque.
    void wrtSineSphericalHarmonicCoefficientsOfBodyUndergoingTorque( Eigen::MatrixXd& partialMatrix,
                                                                     const std::vector< std::pair< int, int > >& blockIndices );

    //! Partial w.r.t. cosine SH coefficient block of body exerting torque.
    void wrtCosineSphericalHarmonicCoefficientsOfBodyExertingTorque( Eigen::MatrixXd& partialMatrix,
                                                                     const std::vector< std::pair< int, int > >& blockIndices );

    //! Partial w.r.t. sine SH coefficient block of body exerting torque.
    void wrtSineSphericalHarmonicCoefficientsOfBodyExertingTorque( Eigen::MatrixXd& partialMatrix,
                                                                   const std::vector< std::pair< int, int > >& blockIndices );

    //! Add body-2 spin torque contribution partial w.r.t. one body-1 coefficient.
    void addBody2SpinTorquePartialWrtBody1Coefficient( Eigen::Vector3d& partial,
                                                       const int degree,
                                                       const int order,
                                                       const bool wrtCosineCoefficient ) const;

    //! Add body-2 spin torque contribution partial w.r.t. one body-2 coefficient.
    void addBody2SpinTorquePartialWrtBody2Coefficient( Eigen::Vector3d& partial,
                                                       const int degree,
                                                       const int order,
                                                       const bool wrtCosineCoefficient );

    //! Partial w.r.t. an effective gravitational parameter used by the associated acceleration model.
    void wrtGravitationalParameter( Eigen::MatrixXd& partialMatrix );

    //! Partial w.r.t. the mass of the body undergoing torque.
    void wrtBodyUndergoingTorqueMass( Eigen::MatrixXd& partialMatrix );

    //! Partial w.r.t. a rotation-model parameter of one body.
    void wrtRotationModelParameter( Eigen::MatrixXd& partialMatrix,
                                    const bool wrtBodyUndergoingTorque,
                                    const estimatable_parameters::EstimatebleParametersEnum parameterType,
                                    const std::string& secondaryIdentifier );

    //! Partial w.r.t. polynomial gravity-field variation amplitudes of one body.
    void wrtPolynomialGravityFieldVariations( const bool wrtBodyUndergoingTorque,
                                              const std::vector< std::pair< int, int > >& cosineBlockIndices,
                                              const std::vector< std::pair< int, int > >& sineBlockIndices,
                                              const std::vector< std::vector< std::pair< int, int > > > powersPerCosineBlockIndex,
                                              const std::vector< std::vector< std::pair< int, int > > > powersPerSineBlockIndex,
                                              const double referenceEpoch,
                                              Eigen::MatrixXd& partialMatrix );

    //! Partial w.r.t. periodic gravity-field variation amplitudes of one body.
    void wrtPeriodicGravityFieldVariations( const bool wrtBodyUndergoingTorque,
                                            const std::vector< std::pair< int, int > >& cosineBlockIndices,
                                            const std::vector< std::pair< int, int > >& sineBlockIndices,
                                            const std::vector< std::vector< std::pair< int, int > > > periodsPerCosineBlockIndex,
                                            const std::vector< std::vector< std::pair< int, int > > > periodsPerSineBlockIndex,
                                            const std::vector< double >& frequencies,
                                            const double referenceEpoch,
                                            Eigen::MatrixXd& partialMatrix );

    std::shared_ptr< gravitation::FullTwoBodySphericalHarmonicTorque > torqueModel_;
    std::shared_ptr< FullTwoBodySphericalHarmonicsGravityPartial > accelerationPartial_;
    std::shared_ptr< gravitation::FullTwoBodySphericalHarmonicAcceleration > accelerationModel_;
    std::shared_ptr< gravitation::EffectiveMutualSphericalHarmonicsField > effectiveMutualPotentialField_;

    std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinationsToUse_;
    std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > body2TorqueCombinationsToUse_;

    Eigen::Matrix< double, 3, 4 > currentPartialWrtQuaternionOfBodyUndergoingTorque_;
    Eigen::Matrix< double, 3, 4 > currentPartialWrtQuaternionOfBodyExertingTorque_;
    Eigen::Matrix3d currentPartialWrtPositionOfBodyUndergoingTorque_;
    Eigen::Matrix3d currentPartialWrtPositionOfBodyExertingTorque_;
    Eigen::Matrix3d currentRotationFromInertialToBodyFixedFrameOfBodyUndergoingTorque_;
    Eigen::Matrix3d currentRotationFromInertialToBodyFixedFrameOfBodyExertingTorque_;
    Eigen::Matrix3d currentRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorque_;
    Eigen::Vector3d currentRelativePositionInInertialFrame_;
    Eigen::Vector3d currentRelativePositionInBodyFixedFrameOfBodyUndergoingTorque_;
    Eigen::Vector3d currentMutualPotentialGradientInBodyFixedFrameOfBodyUndergoingTorque_;
    Eigen::Matrix3d currentBody2SpinTorquePartialWrtBodyFixedRelativePosition_;
    Eigen::MatrixXd currentTransformedCosineCoefficientsBody2_;
    Eigen::MatrixXd currentTransformedSineCoefficientsBody2_;
    std::array< Eigen::MatrixXd, 3 > currentTransformedCosineCoefficientsBody2AngularMomentum_;
    std::array< Eigen::MatrixXd, 3 > currentTransformedSineCoefficientsBody2AngularMomentum_;

    Eigen::MatrixXd body2CoefficientBasisCosineScratch_;
    Eigen::MatrixXd body2CoefficientBasisSineScratch_;
    Eigen::MatrixXd transformedCosineBody2CoefficientPartialsScratch_;
    Eigen::MatrixXd transformedSineBody2CoefficientPartialsScratch_;
    std::array< Eigen::MatrixXd, 3 > transformedCosineCoefficientsBody2AngularMomentumPartialsScratch_;
    std::array< Eigen::MatrixXd, 3 > transformedSineCoefficientsBody2AngularMomentumPartialsScratch_;

    Eigen::MatrixXd partialOfTransformedCosineCoefficientsBody2Scratch_;
    Eigen::MatrixXd partialOfTransformedSineCoefficientsBody2Scratch_;
    std::array< Eigen::MatrixXd, 3 > partialOfTransformedCosineCoefficientsBody2AngularMomentumScratch_;
    std::array< Eigen::MatrixXd, 3 > partialOfTransformedSineCoefficientsBody2AngularMomentumScratch_;
    std::array< std::vector< Eigen::MatrixXcd >, 4 > derivativeOfWignerDMatricesWrtRelativeQuaternionScratch_;

    double currentDistance_;
    double currentCosineOfLatitude_;
    double currentPreMultiplier_;
    double currentBodyUndergoingTorqueMass_;
    std::vector< double > currentRadius1Powers_;
    std::vector< double > currentRadius2Powers_;

    observation_partials::RotationMatrixPartialNamedList rotationMatrixPartialsOfBodyUndergoingTorque_;
    observation_partials::RotationMatrixPartialNamedList rotationMatrixPartialsOfBodyExertingTorque_;
};

}  // namespace acceleration_partials

}  // namespace tudat

#endif  // TUDAT_FULLTWOBODYSPHERICALHARMONICGRAVITATIONALTORQUEPARTIAL_H
