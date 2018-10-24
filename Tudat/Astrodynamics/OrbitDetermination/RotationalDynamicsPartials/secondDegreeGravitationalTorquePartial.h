/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SECONDDEGREEGRAVITATIONALTORQUEPARTIALS_H
#define TUDAT_SECONDDEGREEGRAVITATIONALTORQUEPARTIALS_H

#include "Tudat/Astrodynamics/Gravitation/secondDegreeGravitationalTorque.h"
#include "Tudat/Astrodynamics/OrbitDetermination/RotationalDynamicsPartials/torquePartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/RotationalDynamicsPartials/inertiaTensorPartial.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

namespace tudat
{

namespace acceleration_partials
{

//! Function that computes the partial of degree 2 torque w.r.t. the current quaternion elements
/*!
 * Function that computes the partial of degree 2 torque w.r.t. the current quaternion elements
 * \param premultiplier Model premultiplier
 * \param inertiaTensor Body inertial tensor
 * \param bodyFixedRelativePosition Position of body exerting torque, w.r.t. body undergoing torque,
 * expressed in frame fixed to body undergoing torque
 * \param intertialRelativePosition Position of body exerting torque, w.r.t. body undergoing torque,
 * expressed in inertial frame
 * \param derivativeOfRotationMatrixWrtQuaternions Current partial derivative of rotation matrix from body-fixed to base
 * frame w.r.t. quaternion
 * \return Partial derivative of degree 2 torque w.r.t. the current quaternion elements
 */
Eigen::Matrix< double, 3, 4 > getPartialDerivativeOfSecondDegreeGravitationalTorqueWrtQuaternion(
        const double premultiplier,
        const Eigen::Matrix3d& inertiaTensor,
        const Eigen::Vector3d& bodyFixedRelativePosition,
        const Eigen::Vector3d& intertialRelativePosition,
        const std::vector< Eigen::Matrix3d > derivativeOfRotationMatrixWrtQuaternions );

//! Class to calculate the partials of the central gravitational acceleration w.r.t. parameters and states.
class SecondDegreeGravitationalTorquePartial: public TorquePartial
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param torqueModel Torque model for which partial derivatives are to be computed
     * \param getInertiaTensorNormalizationFactor Function the inertia tensor normalization factor
     * \param acceleratedBody Name of body undergoing torque
     * \param acceleratingBody Name of body exerting torque
     */
    SecondDegreeGravitationalTorquePartial(
            const std::shared_ptr< gravitation::SecondDegreeGravitationalTorqueModel > torqueModel,
            const std::function< double( ) > getInertiaTensorNormalizationFactor,
            const std::string acceleratedBody,
            const std::string acceleratingBody ):
        TorquePartial( acceleratedBody, acceleratingBody, basic_astrodynamics::second_order_gravitational_torque ),
        torqueModel_( torqueModel ), getInertiaTensorNormalizationFactor_( getInertiaTensorNormalizationFactor )
    {
        currentRotationMatrixDerivativesWrtQuaternion_.resize( 4 );
    }

    //! Destructor
    ~SecondDegreeGravitationalTorquePartial( ){ }

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
    void wrtOrientationOfAcceleratedBody(
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

    //! Function to compute partial of torque w.r.t. gravitational parameter
    /*!
     * Function to compute partial of torque w.r.t. gravitational parameter
     * \param gravitationalParameterPartial Computed partial of torque w.r.t. gravitational parameter (returned by reference)
     */
    void wrtGravitationalParameterOfCentralBody( Eigen::MatrixXd& gravitationalParameterPartial );

    //! Function to compute partial of torque w.r.t. spherical harmonic cosine coefficients
    /*!
     * Function to compute partial of torque w.r.t. spherical harmonic cosine coefficients
     * \param sphericalHarmonicCoefficientPartial Computed partial of torque w.r.t. spherical harmonic cosine coefficients
     * (returned by reference)
     * \param c20Index for degree=2,order=0 coefficient
     * \param c21Index for degree=2,order=1 coefficient
     * \param c22Index for degree=2,order=2 coefficient
     */
    void wrtCosineSphericalHarmonicCoefficientsOfCentralBody(
            Eigen::MatrixXd& sphericalHarmonicCoefficientPartial,
            const int c20Index, const int c21Index, const int c22Index );

    //! Function to compute partial of torque w.r.t. spherical harmonic sine coefficients
    /*!
     * Function to compute partial of torque w.r.t. spherical harmonic sine coefficients
     * \param sphericalHarmonicCoefficientPartial Computed partial of torque w.r.t. spherical harmonic sine coefficients
     * (returned by reference)
     * \param s21Index for degree=2,order=1 coefficient
     * \param s22Index for degree=2,order=2 coefficient
     */
    void wrtSineSphericalHarmonicCoefficientsOfCentralBody(
            Eigen::MatrixXd& sphericalHarmonicCoefficientPartial,
            const int s21Index, const int s22Index );

    std::shared_ptr< gravitation::SecondDegreeGravitationalTorqueModel > torqueModel_;

    //! Function the inertia tensor normalization factor
    std::function< double( ) > getInertiaTensorNormalizationFactor_;

    //! Current quaternion elements
    Eigen::Vector4d currentQuaternionVector_;

    //! Current spherical harmonic coefficient pre-multiplier
    Eigen::Matrix3d currentCoefficientPartialPremultiplier_;

    //! Current reltive position of body undergoing acceleration w.r.t. that exerting acceleration, in body-fixed frame
    Eigen::Vector3d currentBodyFixedRelativePosition_;

    //! Current partial derivatives of rotation matrix from body-fixed to base frame w.r.t. quaternion
    std::vector< Eigen::Matrix3d > currentRotationMatrixDerivativesWrtQuaternion_;

    //! Current partial derivative of torque w.r.t. quaternion
    Eigen::Matrix< double, 3, 4 > currentPartialDerivativeWrtQuaternion_;

};

} // namespace acceleration_partials

} // namespace tudat

#endif // TUDAT_SECONDDEGREEGRAVITATIONALTORQUEPARTIALS_H
