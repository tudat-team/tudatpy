/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SPHERICALHARMONICGRAVITATIONALTORQUEPARTIALS_H
#define TUDAT_SPHERICALHARMONICGRAVITATIONALTORQUEPARTIALS_H

#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicGravitationalTorque.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/sphericalHarmonicAccelerationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/RotationalDynamicsPartials/torquePartial.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

namespace tudat
{

namespace acceleration_partials
{

//! Function to compute partial derivative of spherical harmonic torque w.r.t. quaternion elements
/*!
 * Function to compute partial derivative of spherical harmonic torque w.r.t. quaternion elements
 * \param bodyFixedRelativePositionCrossProductMatrix Cross-product matrix for body-fixed relative position of body undergoing
 * torque (w.r.t. body exerting torque)
 * \param bodyFixedPotentialGradientPositionPartial Gradient of spherical harmonic potential, in body-fixed coordinates
 * \param bodyFixedPotentialGradientCrossProductMatrix Cross-product matrix of body-fixed potential gradient
 * \param inertialRelativePosition Inertial relative position of body undergoing
 * torque (w.r.t. body exerting torque)
 * \param derivativeOfRotationMatrixWrtQuaternions Current partial derivative of rotation matrix from body-fixed to base
 * frame w.r.t. quaternion
 * \return Partial derivative of spherical harmonic torque w.r.t. the current quaternion elements
 */
Eigen::Matrix< double, 3, 4 > getPartialDerivativeOfSphericalHarmonicGravitationalTorqueWrtQuaternion(
        const Eigen::Matrix3d& bodyFixedRelativePositionCrossProductMatrix,
        const Eigen::Matrix3d& bodyFixedPotentialGradientPositionPartial,
        const Eigen::Matrix3d& bodyFixedPotentialGradientCrossProductMatrix,
        const Eigen::Vector3d& inertialRelativePosition,
        const std::vector< Eigen::Matrix3d > derivativeOfRotationMatrixWrtQuaternions );

//! Class to calculate the partials of the central gravitational torque w.r.t. parameters and states.
class SphericalHarmonicGravitationalTorquePartial: public TorquePartial
{
public:

    SphericalHarmonicGravitationalTorquePartial(
            const std::shared_ptr< gravitation::SphericalHarmonicGravitationalTorqueModel > torqueModel,
            const std::shared_ptr< acceleration_partials::SphericalHarmonicsGravityPartial > accelerationPartial,
            const std::string acceleratedBody,
            const std::string acceleratingBody ):
        TorquePartial( acceleratedBody, acceleratingBody, basic_astrodynamics::spherical_harmonic_gravitational_torque ),
        torqueModel_( torqueModel ), accelerationPartial_( accelerationPartial )
    {
        currentRotationMatrixDerivativesWrtQuaternion_.resize( 4 );
    }

    ~SphericalHarmonicGravitationalTorquePartial( ){ }

    //! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current torque.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency, 1 otherwise).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
    getParameterPartialFunction( std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter );

    //! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current torque.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter );

    //! Function for calculating the partial of the torque w.r.t. the orientation of the accelerated body.
    /*!
     *  Function for calculating the partial of the torque w.r.t. the orientation of the accelerated body and
     *  adding it to the existing partial block.
     *  \param partialMatrix Block of partial derivatives of torque w.r.t. Orientation of body
     *  undergoing torque where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    virtual void wrtOrientationOfAcceleratedBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 4 ) += currentPartialDerivativeWrtQuaternion_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 4 ) -= currentPartialDerivativeWrtQuaternion_;
        }
    }

    //! Function for determining if the torque is dependent on a non-rotational integrated state.
    /*!
     *  Function for determining if the torque is dependent on a non-rotational integrated state.
     *  \param stateReferencePoint Reference point id of propagated state
     *  \param integratedStateType Type of propagated state for which dependency is to be determined.
     *  \return True if dependency exists (non-zero partial), false otherwise.
     */
    bool isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType )
    {
        bool isStateDerivativeDependent = 0;
        if( ( ( stateReferencePoint.first == bodyUndergoingTorque_ || ( stateReferencePoint.first == bodyExertingTorque_ ) )
              && integratedStateType == propagators::translational_state ) )
        {
            isStateDerivativeDependent = true;
        }
        else if( ( ( stateReferencePoint.first == bodyUndergoingTorque_ || ( stateReferencePoint.first == bodyExertingTorque_ ) )
              && integratedStateType == propagators::body_mass_state ) )
        {
            throw std::runtime_error( "Warning, dependency of 2nd degree gravity torques on body masses not yet implemented" );
        }
        return isStateDerivativeDependent;
    }

    //! Function for calculating the partial of the torque w.r.t. a non-rotational integrated state
    /*!
     *  Function for calculating the partial of the torque w.r.t. a non-rotational integrated state
     *  and adding it to the existing partial block.
     *  \param partialMatrix Block of partial derivatives of where current partial is to be added.
     *  \param stateReferencePoint Reference point id of propagated state
     *  \param integratedStateType Type of propagated state for which partial is to be computed.
     */
    void wrtNonRotationalStateOfAdditionalBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType );

    //! Update partial model to current time
    /*!
     * Update partial model to current time, retrieving information from consituent functions
     * \param currentTime Time to which partial is to be updated
     */
    void update( const double currentTime = TUDAT_NAN );

protected:

    //! Function to compute torque partial from constituent spherical harmonic acceleration partial
    /*!
     * Function to compute torque partial from constituent spherical harmonic acceleration partial
     * \param partialMatrix Matrux with partial derivative (returned by reference)
     * \param accelerationPartialFunction Pair with (first: function returning partial by reference; second: number of columns in
     * partial derivative).
     */
    void getParameterPartialFromAccelerationPartialFunction(
            Eigen::MatrixXd& partialMatrix,
            const std::pair< std::function< void( Eigen::MatrixXd& ) >, int >& accelerationPartialFunction );

    //! Current quaternion elements
    Eigen::Vector4d currentQuaternionVector_;

    //! Current rotation matrix from inertial to body-fixed frame
    Eigen::Matrix3d currentRotationToBodyFixedFrame_;

    //! Current body-fixed relative position of body undergoing torque (w.r.t. body exerting torque)
    Eigen::Vector3d currentBodyFixedRelativePosition_;

    //! Current gradient of spherical harmonic potential, in body-fixed coordinates
    Eigen::Vector3d currentBodyFixedPotentialGradient_;

    //! Cross-product matrix for body-fixed relative position of body undergoing torque (w.r.t. body exerting torque)
    Eigen::Matrix3d currentBodyFixedRelativePositionCrossProductMatrix_;

    //! Cross-product matrix of body-fixed potential gradient
    Eigen::Matrix3d currentBodyFixedPotentialGradientCrossProductMatrix_;

    //! Current matrix by which to pre-multiply acceleration partial to obtain torque partial
    Eigen::Matrix3d currentParameterPartialPreMultiplier_;

    //! Current partial derivative of torque w.r.t. quaternion
    Eigen::Matrix< double, 3, 4 > currentPartialDerivativeWrtQuaternion_;

    //! Current partial derivative of rotation matrix from body-fixed to base frame w.r.t. quaternion
    std::vector< Eigen::Matrix3d > currentRotationMatrixDerivativesWrtQuaternion_;

    //! Torque model for which to compute partials
    std::shared_ptr< gravitation::SphericalHarmonicGravitationalTorqueModel > torqueModel_;

    //! Partial for associated spherical harmonic acceleration
    const std::shared_ptr< acceleration_partials::SphericalHarmonicsGravityPartial > accelerationPartial_;

};

} // namespace acceleration_partials

} // namespace tudat

#endif // TUDAT_SPHERICALHARMONICGRAVITATIONALTORQUEPARTIALS_H
